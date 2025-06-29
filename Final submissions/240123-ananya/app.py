import streamlit as st
import pandas as pd
import numpy as np
import tensorflow as tf
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.io.wavfile import write
from scipy.signal import fftconvolve
from mido import Message, MidiFile, MidiTrack
import pretty_midi
import matplotlib.pyplot as plt
import io


# --------------------------
# CONFIG
# --------------------------
MODEL_PATH = "molecule_dj_model.h5"
CSV_PATH = "name_smiles.csv"
NOTE_MIN, NOTE_MAX = 48, 84
SEQ_LENGTH = 50

# --------------------------
# LOADERS
# --------------------------
@st.cache_resource
def load_model():
    return tf.keras.models.load_model(MODEL_PATH, compile=False)

@st.cache_data
def load_data():
    return pd.read_csv(CSV_PATH)

model = load_model()
df = load_data()

# --------------------------
# UTILS
# --------------------------
def extract_features(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=128)
    return np.array(fp)

def smiles_to_seed(smiles, features):
    melody = []
    base = NOTE_MIN + int(np.sum(features) % (NOTE_MAX - NOTE_MIN - 12))
    scale = [0, 2, 4, 5, 7, 9, 11, 12]
    for i, ch in enumerate(smiles):
        idx = ord(ch) % len(scale)
        note = base + scale[idx]
        note = min(max(NOTE_MIN, note), NOTE_MAX)
        melody.append(note)
    while len(melody) < SEQ_LENGTH:
        melody.append(base + scale[len(melody) % len(scale)])
    return melody[:SEQ_LENGTH]

def generate_music(seed_notes, steps=100):
    output = seed_notes.copy()
    for _ in range(steps):
        X = np.array(output[-SEQ_LENGTH:]).reshape(-1, SEQ_LENGTH, 1) / 127.0
        pred = model.predict(X, verbose=0)
        next_note = np.argmax(pred, axis=1)[0] + NOTE_MIN
        output.append(next_note)
    return output

def add_reverb(signal, sr=44100, decay=0.4):
    ir = np.zeros(int(sr * decay))
    ir[0] = 1.0
    ir[int(0.01 * sr)] = 0.5
    ir[int(0.02 * sr)] = 0.25
    return fftconvolve(signal, ir, mode="full")[:len(signal)]

def save_audio(note_sequence, filename="output.wav"):
    sr = 44100
    total_duration = 36  # ⏰ 36 seconds
    samples_per_note = int(sr * total_duration / len(note_sequence))
    audio = np.array([], dtype=np.float32)

    for note in note_sequence:
        freq = 440.0 * (2 ** ((note - 69) / 12))
        t = np.linspace(0, samples_per_note / sr, samples_per_note, False)
        # Multi-oscillator
        sine = 0.5 * np.sin(2 * np.pi * freq * t)
        saw = 0.3 * (2 * (t * freq - np.floor(0.5 + t * freq)))
        square = 0.2 * np.sign(np.sin(2 * np.pi * freq * t))
        wave = sine + saw + square
        # Fade in/out
        fade_len = int(0.05 * sr)
        envelope = np.ones_like(wave)
        envelope[:fade_len] = np.linspace(0, 1, fade_len)
        envelope[-fade_len:] = np.linspace(1, 0, fade_len)
        wave *= envelope
        audio = np.concatenate([audio, wave])

    audio = add_reverb(audio)
    audio = np.int16(audio / np.max(np.abs(audio)) * 32767)
    write(filename, sr, audio)

def save_midi(note_sequence, filename="output.mid"):
    mid = MidiFile()
    track = MidiTrack()
    mid.tracks.append(track)
    track.append(Message('program_change', program=0, time=0))
    for note in note_sequence:
        track.append(Message('note_on', note=note, velocity=64, time=0))
        track.append(Message('note_off', note=note, velocity=64, time=200))
    mid.save(filename)

def show_piano_roll(midi_file="output.mid"):
    midi_data = pretty_midi.PrettyMIDI(midi_file)
    piano_roll = midi_data.instruments[0].get_piano_roll(fs=100)

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.imshow(piano_roll, aspect='auto', cmap='magma', origin='lower')
    ax.set_xlabel("Time")
    ax.set_ylabel("MIDI Note")
    ax.set_title("MIDI Piano Roll")
    st.pyplot(fig)


# --------------------------
# UI
# --------------------------
st.set_page_config(
    page_title="🎵 Molecule DJ",
    page_icon="🎹",
    layout="centered"
)

st.markdown(
    """
    <style>
    html, body, .stApp {
      background: url('https://static.vecteezy.com/system/resources/previews/021/916/442/non_2x/hexagonal-with-glowing-particles-on-dark-blue-background-science-technology-medicine-chemistry-data-network-background-design-illustration-vector.jpg') no-repeat center center fixed;
      background-size: cover;
      background-color: transparent;
    }

    .block-container {
      background-color: rgba(0, 0, 0, 0); /* transparent to show your bg */
    }

    .css-18e3th9 {
      background-color: rgba(0, 0, 0, 0); /* also for container */
    }

    .audio-container {
        border: 2px solid #00f7ff;
        border-radius: 10px;
        padding: 10px;
        box-shadow: 0 0 20px #00f7ff, 0 0 40px #00f7ff, 0 0 60px #00f7ff;
        background-color: rgba(0,0,0,0.2);
    }
    audio {
        width: 100%;
    }
    </style>
    """,
    unsafe_allow_html=True
)



st.title("🎹 Molecule DJ")
st.markdown(
    "✨ *Generate unique melodies from molecules!*"
)

options = df["Name"].dropna().unique().tolist()
selected_name = st.selectbox("🔬 **Select a molecule**", options)

if st.button("🎼 **Generate Music**"):
    smiles = df[df["Name"] == selected_name]["SMILES"].values[0]
    features = extract_features(smiles)
    if features is not None:
        seed = smiles_to_seed(smiles, features)
        notes = generate_music(seed)
        save_audio(notes)
        save_midi(notes)

        st.success(f"✅ *Music generated for* **{selected_name}**!")

        show_piano_roll("output.mid")  # <- correct indent!

        st.audio("output.wav")

        with open("output.wav", "rb") as f:
            st.download_button(
                "💾 Download WAV",
                f,
                file_name=f"{selected_name}_music.wav",
                mime="audio/wav"
            )

        with open("output.mid", "rb") as f:
            st.download_button(
                "💾 Download MIDI",
                f,
                file_name=f"{selected_name}_music.mid",
                mime="audio/midi"
            )
    else:
        st.error("❌ Invalid SMILES. Please try another molecule.")


st.markdown("---")
st.markdown("*Made with ❤️*")
