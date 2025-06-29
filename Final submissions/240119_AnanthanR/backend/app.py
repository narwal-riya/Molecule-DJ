from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from FinalModel import extract_features, generate_music
from tensorflow.keras.models import load_model
import numpy as np
import os
import logging

# Load model
model = load_model("molecule_dj_model.h5", compile = False)

# Logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = FastAPI()

# CORS setup
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Request/Response schemas
class Prompt(BaseModel):
    smiles: str

class NotesResponse(BaseModel):
    notes: list

@app.post("/generate", response_model=NotesResponse)
def generate_from_smiles(prompt: Prompt):
    features = extract_features(prompt.smiles)
    if features is None:
        return {"notes": []}
    
    # No scaler used here for simplicity — match your model logic
    seed_note = int(60 + (np.sum(features) % 25))
    seed_sequence = [seed_note] * 50
    notes = generate_music(model, seed_sequence)
    return {"notes": notes if notes else []}

@app.get("/health-check")
def health():
    return {"status": "healthy"}

if __name__ == "__main__":
    import uvicorn
    uvicorn.run("app:app", host="0.0.0.0", port=8000, reload=True)
