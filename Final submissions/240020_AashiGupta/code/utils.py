from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
import numpy as np

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
    except Exception:
        return None

def smiles_to_melody(smiles, features, seq_length=SEQ_LENGTH+1):
    scales = [
        [0, 2, 4, 5, 7, 9, 11, 12],      # Major
        [0, 2, 3, 5, 7, 8, 10, 12],      # Minor
        [0, 3, 5, 7, 10, 12],            # Blues
        [0, 2, 5, 7, 9, 12],             # Pentatonic
    ]
    scale_idx = (int(features[1]) + int(features[2])) % len(scales)
    scale = scales[scale_idx]
    key_base = NOTE_MIN + int((np.abs(features[0]) + features[5]) % (NOTE_MAX - NOTE_MIN - 12))
    melody = []
    for i, ch in enumerate(smiles):
        idx = (ord(ch) + int(features[2]) + i*3) % len(scale)
        fp_offset = sum([int(features[6 + ((i+j) % (len(features)-6))]) for j in range(3)]) % 12
        note = key_base + scale[idx] + fp_offset
        note = min(max(NOTE_MIN, note), NOTE_MAX)
        melody.append(note)
    if len(melody) < seq_length:
        for j in range(seq_length - len(melody)):
            melody.append(key_base + scale[j % len(scale)])
    return melody[:seq_length]

def describe_molecule(features):
    return f"Generated from {len(features)} molecular features including mass, polarity, and fingerprint."

def explain_music_mapping(smiles, features, notes):
    scales = [
        "Major (C D E F G A B C)",
        "Minor (C D Eb F G Ab Bb C)",
        "Blues (C Eb F G Bb C)",
        "Pentatonic (C D F G A C)"
    ]
    scale_idx = (int(features[1]) + int(features[2])) % len(scales)
    key_base = NOTE_MIN + int((np.abs(features[0]) + features[5]) % (NOTE_MAX - NOTE_MIN - 12))
    explanation = [
        "🎵 How is this molecule turned into music?",
        "",
        "Step 1: Calculate molecular features",
        f"  • Molecular weight: {features[0]:.2f}",
        f"  • Atom count: {features[1]}",
        f"  • Bond count: {features[2]}",
        f"  • H-bond donors: {features[3]}",
        f"  • H-bond acceptors: {features[4]}",
        f"  • TPSA (surface area): {features[5]:.2f}",
        "",
        "Step 2: Choose a musical key and scale",
        f"  • Scale: {scales[scale_idx]} (chosen by atom and bond count)",
        f"  • Base note: MIDI {key_base} (chosen by molecular weight and TPSA)",
        "",
        "Step 3: Map SMILES characters to notes",
        "  • Each character is mapped using ASCII code, bond count, fingerprint, and position for maximum variety.",
        "",
        "Step 4: Build the melody",
        f"  • First 5 notes: {notes[:5]}",
        "",
        f"SMILES used: {smiles}"
    ]
    return "\n".join(explanation)

def generate_music(model, smiles, features, length=50, temperature=1.0):
    sequence = smiles_to_melody(smiles, features, seq_length=50)
    generated = []
    for _ in range(length):
        input_seq = np.array(sequence[-50:]).reshape(1, 50, 1) / 127.0
        prediction = model.predict(input_seq, verbose=0)[0]
        preds = np.log(prediction + 1e-8) / temperature
        exp_preds = np.exp(preds)
        probs = exp_preds / np.sum(exp_preds)
        next_note_class = np.random.choice(N_NOTES, p=probs)
        next_note = NOTE_MIN + next_note_class
        generated.append(next_note)
        sequence.append(next_note)
    return generated
