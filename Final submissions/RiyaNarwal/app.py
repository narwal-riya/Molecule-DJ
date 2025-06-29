import streamlit as st
import numpy as np
import pretty_midi
import subprocess
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
from tensorflow.keras.models import load_model

model = load_model('molecule_dj_model.keras', compile=False)

st.set_page_config(page_title="Molecule DJ", page_icon="🎵", layout="centered")

st.markdown(
    """
    <style>
    body {background-color: #87CEEB;}
    .main {background-color: #87CEEB;}
    h1, h2, h3 {color: #004466;}
    .stButton > button {
        background-color: #004466;
        color: white;
        border-radius: 8px;
        padding: 10px 20px;
        font-weight: bold;
    }
    </style>
    """,
    unsafe_allow_html=True
)

st.title("🎵 Molecule DJ - SMILES to Music 🎵")
st.subheader("✨ Turn your molecule into a musical sequence!")

st.markdown(
    """
    Enter a **SMILES** string below to generate a *unique* musical sequence.
    This uses a trained **LSTM model** to predict MIDI note patterns based on molecular features.
    """,
    unsafe_allow_html=True
)

def extract_features(smiles):
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

def generate_music(seed_sequence, length=50):
    generated = seed_sequence.copy()
    for _ in range(length):
        X_pred = np.array(generated[-10:]).reshape(1, 10, 1)
        next_note = model.predict(X_pred, verbose=0)[0][0]
        generated.append(int(next_note))
    return generated

def notes_to_midi(note_sequence, filename="output.mid"):
    midi = pretty_midi.PrettyMIDI()
    piano = pretty_midi.Instrument(program=0)
    time = 0
    for note_num in note_sequence:
        note = pretty_midi.Note(velocity=100, pitch=int(note_num), start=time, end=time+0.5)
        piano.notes.append(note)
        time += 0.5
    midi.instruments.append(piano)
    midi.write(filename)
    return filename

def midi_to_wav(midi_path, soundfont_path, wav_path):
    subprocess.run([
        "fluidsynth", "-ni", soundfont_path, midi_path, "-F", wav_path, "-r", "44100"
    ])
    return wav_path

st.markdown("---")
smiles_input = st.text_input("🧪 Enter your molecule's SMILES string:")

if st.button("🎶 Generate Music"):
    if not smiles_input:
        st.warning("⚠️ Please enter a SMILES string to continue!")
    else:
        features = extract_features(smiles_input)
        if features is None:
            st.error("❌ Invalid SMILES string. Please try again!")
        else:
            features_scaled = (features - features.mean()) / features.std()
            note = int(60 + (np.sum(features_scaled) % 25))
            seed_sequence = [note] * 10
            music_sequence = generate_music(seed_sequence)

            st.success("✅ Generated Musical Sequence!")
            st.code(music_sequence, language='python')

            midi_file = "generated.mid"
            wav_file = "generated.wav"
            soundfont_file = "soundfont.sf2"

            notes_to_midi(music_sequence, midi_file)
            midi_to_wav(midi_file, soundfont_file, wav_file)

            with open(wav_file, "rb") as f:
                audio_bytes = f.read()
            st.audio(audio_bytes, format='audio/wav')
