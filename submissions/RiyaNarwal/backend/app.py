from fastapi import FastAPI
from fastapi.responses import StreamingResponse, JSONResponse
from pydantic import BaseModel
from fastapi.middleware.cors import CORSMiddleware
from gtts import gTTS
from io import BytesIO
from music_generator import generate_audio_from_text
import uvicorn

# Create FastAPI app
app = FastAPI()

# Let frontend use backend
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # Change to ["http://localhost:5173"] in production
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Request and response models
class InputText(BaseModel):
    text: str

class HealthCheckResponse(BaseModel):
    status: str

# Health check (POST)
@app.post("/health-check")
def health_check():
    return HealthCheckResponse(status="healthy")

# Create an audio file from text (POST)
@app.post("/generate")
def generate_audio(req: InputText):
    tts = gTTS(req.text)
    audio_fp = BytesIO()
    tts.write_to_fp(audio_fp)
    audio_fp.seek(0)
    return StreamingResponse(audio_fp, media_type="audio/wav")

# Run app
if __name__ == "__main__":
    uvicorn.run(app, host="127.0.0.4", port=12345)

