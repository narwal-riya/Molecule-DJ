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

# -----------------------------
# 📚 Setup
# -----------------------------

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

NOTE_RANGE = (48, 84)  # MIDI note range: C3 to C6
SEQUENCE_LEN = 50
FP_SIZE = 128  # Morgan fingerprint size

fp_gen = GetMorganGenerator(radius=2, fpSize=FP_SIZE)

# -----------------------------
# 🔬 Functions
# -----------------------------

def get_mol_features(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    # Basic descriptors
    desc = [
        Descriptors.MolWt(mol),
        mol.GetNumAtoms(),
        mol.GetNumBonds(),
        Descriptors.NumHDonors(mol),
        Descriptors.NumHAcceptors(mol),
        Descriptors.TPSA(mol),
    ]
    # Morgan fingerprint bits
    fp = fp_gen.GetFingerprint(mol).ToList()
    return np.array(desc + fp)

def smiles_to_notes(smiles, features, length=SEQUENCE_LEN+1):
    """Turn SMILES into a pseudo-melody sequence."""
    scales = [
        [0, 2, 4, 5, 7, 9, 11, 12],   # Major
        [0, 2, 3, 5, 7, 8, 10, 12],   # Minor
        [0, 3, 5, 6, 7, 10, 12],      # Blues-ish
    ]
    scale_idx = int(features[1] + features[2]) % len(scales)
    base = NOTE_RANGE[0] + int(features[0] % (NOTE_RANGE[1] - NOTE_RANGE[0] - 12))
    notes = []
    scale = scales[scale_idx]

    for i, char in enumerate(smiles):
        idx = (ord(char) + i) % len(scale)
        bit_offset = int(features[6 + (i % (len(features) - 6))]) % 12
        note = base + scale[idx] + bit_offset
        note = max(NOTE_RANGE[0], min(note, NOTE_RANGE[1]))
        notes.append(note)

    while len(notes) < length:
        notes.append(base + scale[len(notes) % len(scale)])

    return notes[:length]

def load_smiles(csv="name_smiles.csv", limit=5000):
    df = pd.read_csv(csv)
    if "SMILES" not in df.columns:
        raise ValueError("CSV must have a 'SMILES' column.")
    return df["SMILES"].dropna().sample(n=min(len(df), limit)).tolist()

def prep_data(smiles_list):
    X, y = [], []
    feats = []
    for smiles in smiles_list:
        f = get_mol_features(smiles)
        if f is None:
            continue
        notes = smiles_to_notes(smiles, f)
        for i in range(len(notes) - SEQUENCE_LEN):
            seq = notes[i : i + SEQUENCE_LEN]
            target = notes[i + SEQUENCE_LEN] - NOTE_RANGE[0]
            X.append(seq)
            y.append(target)
        feats.append(f)

    X = np.array(X).reshape(-1, SEQUENCE_LEN, 1) / 127.0
    y = np.array(y)
    scaler = StandardScaler().fit(feats)
    return X, y, scaler

def build_network():
    model = Sequential([
        LSTM(128, return_sequences=True, input_shape=(SEQUENCE_LEN, 1)),
        Dropout(0.2),
        LSTM(64),
        Dense(64, activation='relu'),
        Dense(NOTE_RANGE[1] - NOTE_RANGE[0] + 1, activation='softmax')
    ])
    model.compile(optimizer='adam', loss='sparse_categorical_crossentropy', metrics=['accuracy'])
    return model

# -----------------------------
# 🚀 Train
# -----------------------------

def train():
    log.info("📥 Reading SMILES...")
    smiles_list = load_smiles(limit=3000)
    log.info(f"✅ Loaded {len(smiles_list)} molecules.")

    log.info("⚙️ Generating sequences...")
    X, y, scaler = prep_data(smiles_list)
    log.info(f"🔢 Training on {X.shape[0]} sequences.")

    model = build_network()
    log.info("🎵 Starting training...")
    model.fit(X, y, epochs=20, batch_size=64, verbose=1)

    log.info("💾 Saving artifacts...")
    model.save("molecule_dj_custom.keras")
    with open("scaler_custom.pkl", "wb") as f:
        pickle.dump(scaler, f)
    log.info("✅ Done! Model and scaler saved.")

# -----------------------------
# ▶️ Run
# -----------------------------

if __name__ == "__main__":
    train()
