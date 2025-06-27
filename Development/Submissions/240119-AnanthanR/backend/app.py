from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import StreamingResponse
from pydantic import BaseModel
from music_generator import generate_audio_from_text
import uvicorn

app = FastAPI()

# Allow CORS for frontend
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

class Prompt(BaseModel):
    prompt: str

class HealthResponse(BaseModel):
    status: str

@app.post("/health-check")
def health_check():
    return HealthResponse(status="healthy")

@app.post("/generate")
def generate_audio(prompt: Prompt):
    audio = generate_audio_from_text(prompt.prompt)
    return StreamingResponse(audio, media_type="audio/wav")

if __name__ == "__main__":
    uvicorn.run("app:app", host="0.0.0.0", port=8000, reload=True)
