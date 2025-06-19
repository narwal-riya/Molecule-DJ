# piano_engine.py
import numpy as np
from scipy.io.wavfile import write
import io
import hashlib
import random

RATE = 44100

NOTE_FREQS = {
    'C': 261.63, 'D': 293.66, 'E': 329.63,
    'F': 349.23, 'G': 392.00, 'A': 440.00, 'B': 493.88
}

CHORDS = {
    "C": ["C", "E", "G"],
    "G": ["G", "B", "D"],
    "Am": ["A", "C", "E"],
    "F": ["F", "A", "C"],
    "Dm": ["D", "F", "A"],
    "Em": ["E", "G", "B"]
}

CHORD_LIST = list(CHORDS.keys())

def envelope(t, duration):
    attack = 0.02
    release = 0.2
    env = np.ones_like(t)
    env[t < attack] = t[t < attack] / attack
    env[t > duration - release] *= np.maximum(0, (duration - t[t > duration - release]) / release)
    return env

def piano_tone(freq, duration, rate=RATE):
    t = np.linspace(0, duration, int(rate * duration), endpoint=False)
    env = envelope(t, duration)

    sine = np.sin(2 * np.pi * freq * t)
    overtone1 = 0.4 * np.sin(2 * np.pi * freq * 2 * t)
    overtone2 = 0.2 * np.sin(2 * np.pi * freq * 3 * t)

    waveform = (sine + overtone1 + overtone2) * env
    return 0.5 * waveform

def chord_wave(notes, duration):
    waves = [piano_tone(NOTE_FREQS[n], duration) for n in notes]
    return sum(waves) / len(waves)

def melody_note_from_char(char, scale):
    index = ord(char.lower()) % len(scale)
    return scale[index]

def generate_musical_melody(text, target_duration=7.0):
    melody = []
    total_duration = 0.0
    i = 0

    # Seed randomness from text
    seed = int(hashlib.sha256(text.encode()).hexdigest(), 16) % (10 ** 8)
    rng = random.Random(seed)

    # Derive scale and chord progression from hash
    scale = ['C', 'D', 'E', 'F', 'G', 'A', 'B']
    rng.shuffle(scale)

    chord_progression = rng.sample(CHORD_LIST, 4)

    while total_duration < target_duration:
        char = text[i % len(text)]
        chord_name = chord_progression[i % len(chord_progression)]
        chord = CHORDS[chord_name]

        # Vary duration slightly
        duration = rng.choice([0.3, 0.4, 0.5])
        total_duration += duration + 0.03

        # Harmony
        harmony = chord_wave(chord, duration)

        # Melody
        melody_note = melody_note_from_char(char, scale)
        melody_wave = piano_tone(NOTE_FREQS[melody_note], duration)

        final = 0.6 * harmony + 0.4 * melody_wave
        melody.append(final)

        rest = np.zeros(int(0.03 * RATE))
        melody.append(rest)
        i += 1

    return np.concatenate(melody)

def generate_audio_from_text(text):
    wave = generate_musical_melody(text, target_duration=7.0)
    wave = (wave * 32767).astype(np.int16)
    buffer = io.BytesIO()
    write(buffer, RATE, wave)
    buffer.seek(0)
    return buffer
