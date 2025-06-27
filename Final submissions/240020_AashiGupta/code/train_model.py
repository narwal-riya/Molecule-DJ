import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import LSTM, Dense, Dropout
from sklearn.preprocessing import StandardScaler
import pickle
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

generator = GetMorganGenerator(radius=2, fpSize=128)
NOTE_MIN, NOTE_MAX = 48, 84  # C3 to C6, musically useful piano range
N_NOTES = NOTE_MAX - NOTE_MIN + 1
SEQ_LENGTH = 50

def extract_features(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        desc = [
            Descriptors.MolWt(mol),
            mol.GetNumAtoms(),
            mol.GetNumBonds(),
            Descriptors.NumHDonors(mol),
            Descriptors.NumHAcceptors(mol),
            Descriptors.TPSA(mol),
        ]
        fp = generator.GetFingerprint(mol).ToList()
        return np.concatenate([desc, fp])
    except Exception as e:
        logger.warning(f"Failed for SMILES {smiles}: {e}")
        return None

def smiles_to_melody(smiles, features, seq_length=SEQ_LENGTH+1):
    scales = [
        [0, 2, 4, 5, 7, 9, 11, 12],      # Major
        [0, 2, 3, 5, 7, 8, 10, 12],      # Minor
        [0, 3, 5, 7, 10, 12],            # Blues
        [0, 2, 5, 7, 9, 12],             # Pentatonic
    ]
    # Use atom count + bond count for scale selection
    scale_idx = (int(features[1]) + int(features[2])) % len(scales)
    scale = scales[scale_idx]
    # Use molecular weight + TPSA for base note, always in 48–84
    key_base = NOTE_MIN + int((np.abs(features[0]) + features[5]) % (NOTE_MAX - NOTE_MIN - 12))
    melody = []
    for i, ch in enumerate(smiles):
        idx = (ord(ch) + int(features[2]) + i*3) % len(scale)
        # Use several fingerprint bits for offset to maximize diversity
        fp_offset = sum([int(features[6 + ((i+j) % (len(features)-6))]) for j in range(3)]) % 12
        note = key_base + scale[idx] + fp_offset
        note = min(max(NOTE_MIN, note), NOTE_MAX)
        melody.append(note)
    if len(melody) < seq_length:
        for j in range(seq_length - len(melody)):
            melody.append(key_base + scale[j % len(scale)])
    return melody[:seq_length]

def load_smiles(csv_path="molecules.csv", max_rows=None):
    df = pd.read_csv(csv_path)
    if "SMILES" not in df.columns:
        raise ValueError("❌ 'SMILES' column not found.")
    return df["SMILES"].dropna().sample(n=min(max_rows, len(df))).tolist()

def preprocess_data(smiles_list, seq_length=SEQ_LENGTH):
    features = []
    valid_smiles = []
    X, y = [], []
    for smiles in smiles_list:
        feat = extract_features(smiles)
        if feat is not None:
            seed = smiles_to_melody(smiles, feat, seq_length + 1)
            for i in range(len(seed) - seq_length):
                X.append(seed[i:i + seq_length])
                y.append(seed[i + seq_length] - NOTE_MIN)
            features.append(feat)
            valid_smiles.append(smiles)
    X = np.array(X).reshape(-1, seq_length, 1) / 127.0
    y = np.array(y)
    scaler = StandardScaler().fit(features)
    return X, y, scaler, valid_smiles

def build_model(seq_length=SEQ_LENGTH):
    model = Sequential([
        LSTM(128, return_sequences=True, input_shape=(seq_length, 1)),
        Dropout(0.3),
        LSTM(64),
        Dense(64, activation='relu'),
        Dense(N_NOTES, activation='softmax')
    ])
    model.compile(optimizer='adam', loss='sparse_categorical_crossentropy', metrics=['accuracy'])
    return model

def train_model():
    logger.info("🚀 Loading SMILES data...")
    smiles_list = load_smiles("molecules.csv", max_rows=16000)
    logger.info(f"🔬 Extracting features and preparing sequences for {len(smiles_list)} molecules...")
    X, y, scaler, valid_smiles = preprocess_data(smiles_list)
    logger.info("🎼 Building model...")
    model = build_model(seq_length=SEQ_LENGTH)
    logger.info("🎹 Training model (this may take a while on CPU)...")
    model.fit(X, y, epochs=20, batch_size=64, verbose=1)
    logger.info("💾 Saving model and scaler...")
    model.save("molecule_dj_model_full.keras")
    with open("scaler_full.pkl", "wb") as f:
        pickle.dump(scaler, f)
    logger.info("✅ Training complete. Model and scaler saved.")

if __name__ == "__main__":
    train_model()
