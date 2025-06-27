from fastapi import FastAPI
from fastapi.responses import StreamingResponse
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from music_generator import generate_audio_from_text
import uvicorn
import io

app = FastAPI()

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

class InputText(BaseModel):
    text: str

class HealthResponse(BaseModel):
    message: str

@app.post("/health-check")
async def health_check():
    return {"message": "healthy"}

@app.post("/generate")
async def generate_audio(req: InputText):
    audio_buffer = generate_audio_from_text(req.text)
    return StreamingResponse(
        io.BytesIO(audio_buffer.read()),
        media_type="audio/wav"
    )

if __name__ == "__main__":
    # Using host 127.0.0.{length_of_name} => 'RiyaNarwal' is 9 letters
    uvicorn.run("app:app", host="127.0.0.9", port=8000, reload=True)
