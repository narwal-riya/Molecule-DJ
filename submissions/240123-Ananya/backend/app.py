from fastapi import FastAPI, Response
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from music_generator import generate_audio_from_text
from starlette.responses import StreamingResponse
import io

app = FastAPI()

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"]
)

class TextRequest(BaseModel):
    text: str

@app.get("/health-check")
def health_check():
    return Response(content="healthy", media_type="text/plain")

@app.post("/generate")
def generate(request: TextRequest):
    audio_data = generate_audio_from_text(request.text)
    return StreamingResponse(audio_data, media_type="audio/wav")

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="127.0.0.1", port=8000)
