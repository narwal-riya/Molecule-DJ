import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, LSTM, Dropout
from sklearn.preprocessing import StandardScaler
import logging

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def load_smiles_data(csv_path="Final submissions/240123-ananya/name_smiles.csv"):
    """Load SMILES data from your custom CSV."""
    try:
        df = pd.read_csv(csv_path)
        smiles_list = df['SMILES'].dropna().tolist()
        logger.info(f"Loaded {len(smiles_list)} SMILES strings from {csv_path}")
        return smiles_list
    except Exception as e:
        logger.error(f"Error loading dataset: {e}")
        raise


def extract_features(smiles):
    """Extract molecular features from a SMILES string."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        features = [
            Descriptors.MolWt(mol),
            mol.GetNumAtoms(),
            mol.GetNumBonds(),
            Descriptors.NumHDonors(mol),
            Descriptors.NumHAcceptors(mol),
            Descriptors.TPSA(mol)
        ]
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=128)
        fp_array = np.array(fp)
        return np.concatenate([features, fp_array])
    except Exception as e:
        logger.warning(f"Invalid SMILES {smiles}: {e}")
        return None


def preprocess_data(smiles_list, seq_length=50):
    """Preprocess SMILES data into features and pseudo-musical sequences."""
    scaler = StandardScaler()
    features_list = []
    valid_smiles = []

    for smiles in smiles_list[:1000]:
        feat = extract_features(smiles)
        if feat is not None:
            features_list.append(feat)
            valid_smiles.append(smiles)

    if not features_list:
        raise ValueError("No valid SMILES strings processed")

    features_array = scaler.fit_transform(features_list)

    musical_sequences = []
    for feat in features_array:
        note = int(60 + (np.sum(feat) % 25))
        sequence = [note] * (seq_length + 10)  # make it longer!
        musical_sequences.append(sequence)

    X = []
    y = []
    for seq in musical_sequences:
        for i in range(0, len(seq) - seq_length):
            X.append(seq[i:i + seq_length])
            y.append(seq[i + seq_length])

    X = np.array(X)
    y = np.array(y)
    X = X.reshape((X.shape[0], X.shape[1], 1))

    logger.info(f"Prepared {len(X)} sequences for training")
    return X, y, scaler, valid_smiles


def build_model(seq_length):
    """Build and compile the LSTM model."""
    model = Sequential([
        LSTM(128, input_shape=(seq_length, 1), return_sequences=True),
        Dropout(0.2),
        LSTM(64),
        Dropout(0.2),
        Dense(32, activation='relu'),
        Dense(1, activation='linear')
    ])
    model.compile(optimizer='adam', loss='mse')
    logger.info("Model built and compiled")
    return model


def train_model(model, X, y, epochs=10, batch_size=32):
    """Train the model."""
    try:
        model.fit(X, y, epochs=epochs, batch_size=batch_size, verbose=1)
        model.save('molecule_dj_model.h5')
        logger.info("Model trained and saved as molecule_dj_model.h5")
    except Exception as e:
        logger.error(f"Error during training: {e}")
        raise


def generate_music(model, seed_sequence, length=100):
    """Generate a musical sequence from a seed."""
    generated = seed_sequence.copy()
    for _ in range(length):
        X_pred = np.array(generated[-50:]).reshape(1, 50, 1)
        next_note = model.predict(X_pred, verbose=0)[0][0]
        generated.append(int(next_note))
    return generated


def main():
    """Main entry point."""
    try:
        smiles_list = load_smiles_data()
        X, y, scaler, valid_smiles = preprocess_data(smiles_list)

        model = build_model(seq_length=50)
        train_model(model, X, y)

        seed_smiles = valid_smiles[0]
        seed_features = extract_features(seed_smiles)
        if seed_features is not None:
            seed_features = scaler.transform([seed_features])[0]
            seed_note = int(60 + (np.sum(seed_features) % 25))
            seed_sequence = [seed_note] * 50
            music_sequence = generate_music(model, seed_sequence)
            logger.info(f"Generated music sequence: {music_sequence[:10]}...")
        else:
            logger.error("Failed to generate features for seed SMILES")

    except Exception as e:
        logger.error(f"Error in main: {e}")
        raise


if __name__ == "__main__":
    main()
