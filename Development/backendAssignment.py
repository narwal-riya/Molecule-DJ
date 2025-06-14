from fastapi import FastAPI
from fastapi.responses import StreamingResponse
from pydantic import BaseModel
from music_generator import generate_audio_from_text
import uvicorn
from fastapi.middleware.cors import CORSMiddleware


# create a server instance of FastAPI


# create classes as template for your request and response like : Class request(BaseModel) ....


# (optional) add a middleware to handle the cross-origin access and allow only your frontend to access the backend


# create the route "health-check" which will return text : "healthy"


# create a route "generate" which takes req as InputText class defined earlier and use the generate_audio_from_text function to generate audio
# and return it as response to frontend


# Run the app on port 8000 with host = "127.0.0.{length_of_your_name}"
