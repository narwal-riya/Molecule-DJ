from flask import Flask, request, jsonify, send_from_directory, send_file
from flask_cors import CORS
import os, tempfile, pickle
import numpy as np
import tensorflow as tf
from midiutil import MIDIFile
from utils import extract_features, describe_molecule, generate_music, explain_music_mapping

app = Flask(__name__, static_folder="frontend")
CORS(app)

model = tf.keras.models.load_model("molecule_dj_model_full.keras")
with open("scaler_full.pkl", "rb") as f:
    scaler = pickle.load(f)

def notes_to_midi(notes, path):
    midi = MIDIFile(1)
    midi.addTempo(0, 0, 120)
    for i, note in enumerate(notes):
        midi.addNote(0, 0, note, i * 0.5, 1, 100)
    with open(path, "wb") as f:
        midi.writeFile(f)

@app.route("/")
def index():
    return send_from_directory("frontend", "index.html")

@app.route("/generate", methods=["POST"])
def generate():
    data = request.get_json()
    smiles = data.get("smiles", "")
    features = extract_features(smiles)

    if features is None:
        return jsonify({"error": "Invalid SMILES"}), 400

    notes = generate_music(model, smiles, features)
    explanation = explain_music_mapping(smiles, features, notes)

    temp_file = tempfile.NamedTemporaryFile(delete=False, suffix=".mid", dir=".")
    notes_to_midi(notes, temp_file.name)
    temp_file.close()

    return jsonify({
        "notes": notes,
        "explanation": explanation,
        "midi_url": f"/get-midi/{os.path.basename(temp_file.name)}"
    })

@app.route("/get-midi/<filename>")
def get_midi(filename):
    path = os.path.join(".", filename)
    return send_file(path, as_attachment=True, download_name="generated.mid")

@app.route("/<path:filename>")
def serve_static(filename):
    return send_from_directory("frontend", filename)

if __name__ == "__main__":
    app.run(debug=True, port=5000)
