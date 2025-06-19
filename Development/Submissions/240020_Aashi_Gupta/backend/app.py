from fastapi import FastAPI
from fastapi.responses import StreamingResponse
from pydantic import BaseModel
from music_generator import generate_audio_from_text
import uvicorn
from fastapi.middleware.cors import CORSMiddleware

app = FastAPI()

app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173"], 
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


class InputText(BaseModel):
    text: str

@app.get("/health-check")
def health_check():
    return "healthy"

@app.post("/generate")
def generate(req: InputText):
    audio_stream = generate_audio_from_text(req.text)
    return StreamingResponse(audio_stream, media_type="audio/wav")

if __name__ == "__main__":
    uvicorn.run(app, host="127.0.0.5", port=8000)




# from fastapi import FastAPI
# from fastapi.responses import StreamingResponse
# from pydantic import BaseModel
# from music_generator import generate_audio_from_text
# import uvicorn
# from fastapi.middleware.cors import CORSMiddleware


# # create a server instance of FastAPI
# app = FastAPI()

# # (optional) add a middleware to handle the cross-origin access and allow only your frontend to access the backend
# app.add_middleware(
#     CORSMiddleware,
#     allow_origins=["http://localhost:3000"],  # frontend URL
#     allow_credentials=True,
#     allow_methods=["*"],
#     allow_headers=["*"],
# )

# # create classes as template for your request and response like : Class request(BaseModel) ....
# class ChatRequest(BaseModel):
#     user_id: str
#     message: str
#     conversation_id: str = None

# # create the route "health-check" which will return text : "healthy"
# @app.get("/health-check")
# def health_check():
#     return "healthy"

# # create a route "generate" which takes req as InputText class defined earlier and use the generate_audio_from_text function to generate audio
# # and return it as response to frontend
# @app.post("/generate")
# def generate(req: ChatRequest):
#     audio_stream = generate_audio_from_text(req.message)
#     return StreamingResponse(audio_stream, media_type="audio/wav")

# # Run the app on port 8000 with host = "127.0.0.{length_of_your_name}"
# if __name__ == "__main__":
#     uvicorn.run(app, host="127.0.0.5", port=8000)




# backend- path, 
# uvicorn app:app --reload --host 127.0.0.5 --port 8000

#frontend- path,
# npm run dev -- --open
