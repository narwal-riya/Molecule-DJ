import kagglehub
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, LSTM, Dropout
from sklearn.preprocessing import StandardScaler
import logging
import os

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def download_dataset():
    """Download the Big Molecules SMILES Dataset from Kaggle."""
    try:
        path = kagglehub.dataset_download("yanmaksi/big-molecules-smiles-dataset")
        logger.info(f"Dataset downloaded to {path}")
        return path
    except Exception as e:
        logger.error(f"Failed to download dataset: {e}")
        raise

def load_smiles_data(dataset_path):
    """Load SMILES data from the dataset."""
    try:
        csv_path = os.path.join(dataset_path, "SMILES_Big_Data_Set.csv")  # Adjust based on actual file name
        df = pd.read_csv(csv_path)
        smiles_list = df['SMILES'].dropna().tolist()  # Adjust column name as needed
        logger.info(f"Loaded {len(smiles_list)} SMILES strings")
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
            Descriptors.MolWt(mol),  # Molecular weight
            mol.GetNumAtoms(),       # Number of atoms
            mol.GetNumBonds(),       # Number of bonds
            Descriptors.NumHDonors(mol),  # Hydrogen bond donors
            Descriptors.NumHAcceptors(mol),  # Hydrogen bond acceptors
            Descriptors.TPSA(mol),   # Topological polar surface area
        ]
        # Add Morgan fingerprint (bit vector)
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

    # Extract features
    for smiles in smiles_list[:1000]:  # Limit for demo purposes
        feat = extract_features(smiles)
        if feat is not None:
            features_list.append(feat)
            valid_smiles.append(smiles)

    if not features_list:
        raise ValueError("No valid SMILES strings processed")

    # Normalize features
    features_array = scaler.fit_transform(features_list)

    # Generate pseudo-musical sequences (MIDI notes 60-84 for simplicity)
    # Map feature magnitude to note range
    musical_sequences = []
    for feat in features_array:
        # Simplified mapping: sum of normalized features to MIDI note
        note = int(60 + (np.sum(feat) % 25))  # MIDI notes 60-84
        # Create a sequence of length seq_length * 2 for sliding window
        sequence = [note] * (seq_length * 2) # Increase sequence length
        musical_sequences.append(sequence)

    # Prepare data for LSTM using a sliding window
    X = []
    y = []
    for seq in musical_sequences:
        # Iterate through the sequence to create input-output pairs using a sliding window
        for i in range(len(seq) - seq_length):
            X.append(seq[i:i + seq_length])
            y.append(seq[i + seq_length])


    X = np.array(X)
    y = np.array(y)
    # Ensure X is not empty before reshaping
    if X.size > 0:
        X = X.reshape((X.shape[0], X.shape[1], 1))  # Reshape for LSTM
    else:
         raise ValueError("No sequences generated for training. Check data and seq_length.")


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
        Dense(1, activation='linear')  # Predict next MIDI note
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
        # Ensure the seed sequence is at least seq_length long
        if len(generated) < 50:
            logger.error("Seed sequence is too short for generation.")
            return None
        X_pred = np.array(generated[-50:]).reshape(1, 50, 1)
        next_note = model.predict(X_pred, verbose=0)[0][0]
        generated.append(int(next_note))
    return generated

def main():
    """Main function to process dataset and train model."""
    try:
        # Download and load dataset
        dataset_path = download_dataset()
        smiles_list = load_smiles_data(dataset_path)

        # Preprocess data
        X, y, scaler, valid_smiles = preprocess_data(smiles_list)

        # Build and train model
        model = build_model(seq_length=50)
        train_model(model, X, y)

        # Example: Generate music for a random molecule
        if valid_smiles:
            seed_smiles = valid_smiles[0]
            seed_features = extract_features(seed_smiles)
            if seed_features is not None:
                seed_features = scaler.transform([seed_features])[0]
                seed_note = int(60 + (np.sum(seed_features) % 25))
                seed_sequence = [seed_note] * 50  # Initial seed sequence of length seq_length
                music_sequence = generate_music(model, seed_sequence)
                if music_sequence:
                    logger.info(f"Generated music sequence: {music_sequence[:10]}...")  # Log first 10 notes
                    return music_sequence
                else:
                     logger.error("Music sequence generation failed.")
                     return None
            else:
                logger.error("Failed to generate features for seed SMILES")
                return None
        else:
             logger.error("No valid smiles found to generate music for.")
             return None
    except Exception as e:
        logger.error(f"Error in main: {e}")
        raise


if __name__ == "__main__":
    main()