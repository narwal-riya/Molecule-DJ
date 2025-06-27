let generatedNotes = [];
let djSynth, part;
let isPlaying = false, progressTimer;
let musicEndTime = 0;

window.addEventListener("DOMContentLoaded", () => {
    document.getElementById("loading").classList.add("hidden");
});

// Dropdown: scrollable if long
document.getElementById("moleculeSelect").addEventListener("change", function() {
    const smilesInput = document.getElementById("smiles");
    if (this.value) smilesInput.value = this.value;
});

async function generateMusic() {
    const smiles = document.getElementById("smiles").value.trim();
    const error = document.getElementById("error");
    const output = document.getElementById("output");
    const history = document.getElementById("history");
    const controls = document.getElementById("controls");
    const downloadLink = document.getElementById("downloadLink");
    const loading = document.getElementById("loading");
    const playPauseBtn = document.getElementById("playPauseBtn");
    const progress = document.getElementById("progressBar");

    error.innerText = "";
    output.innerText = "";
    controls.classList.add("hidden");
    loading.classList.remove("hidden");

    if (!smiles) {
        loading.classList.add("hidden");
        error.innerText = "⚠️ Please enter a valid SMILES string.";
        return;
    }

    try {
        const res = await fetch("/generate", {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify({ smiles })
        });
        const data = await res.json();

        loading.classList.add("hidden");

        if (!res.ok) {
            error.innerText = `❌ ${data.error}`;
            return;
        }

        generatedNotes = data.notes;
        output.innerText = data.explanation + "\n\nNotes: " + generatedNotes.join(", ");
        controls.classList.remove("hidden");
        downloadLink.href = data.midi_url;
        console.log("Generated notes:", generatedNotes);

        const li = document.createElement("li");
        li.textContent = `${smiles} → ${generatedNotes.slice(0, 10).join(", ")}...`;
        history.prepend(li);

        updateDJRoll(generatedNotes);
        setupDJPlayback();

        playPauseBtn.textContent = "▶️ Play DJ";
        progress.value = 0;
    } catch (err) {
        loading.classList.add("hidden");
        error.innerText = "❌ Failed to connect to backend.";
    }
}

function setupDJPlayback() {
    if (djSynth) djSynth.dispose();

    djSynth = new Tone.PolySynth(Tone.Synth, {
        oscillator: { type: "triangle" },
        envelope: {
            attack: 0.01,
            decay: 0.09,
            sustain: 0.15,
            release: 0.13
        }
    });

    const reverb = new Tone.Reverb({ decay: 1.1, wet: 0.09 }).toDestination();
    djSynth.connect(reverb);

    Tone.Transport.stop();
    Tone.Transport.cancel();

    if (part) part.dispose();

    const noteInterval = 0.22;
    const noteDuration = Tone.Time("8n").toSeconds();
    const events = generatedNotes.map((note, i) => {
        let midiNote = note
        return {
            time: i * noteInterval,
            note: Tone.Frequency(midiNote, "midi").toNote(),
            duration: "8n",
            velocity: 0.75
        };
    });
    const lastEvent = events[events.length - 1];
    musicEndTime = (lastEvent ? lastEvent.time + noteDuration : 0);

    part = new Tone.Part((time, value) => {
        djSynth.triggerAttackRelease(value.note, value.duration, time, value.velocity);
    }, events);

    part.loop = false;
    part.start(0);

    document.getElementById("progressBar").value = 0;
    isPlaying = false;
}

async function togglePlayPause() {
    if (Tone.context.state !== "running") {
        await Tone.start();
    }
    const playPauseBtn = document.getElementById("playPauseBtn");
    const progress = document.getElementById("progressBar");
    const duration = musicEndTime;

    if (isPlaying) {
        Tone.Transport.pause();
        playPauseBtn.textContent = "▶️ Play DJ";
        clearInterval(progressTimer);
        isPlaying = false;
    } else {
        Tone.Transport.start();
        playPauseBtn.textContent = "⏸️ Pause DJ";
        progress.max = duration;
        isPlaying = true;

        progressTimer = setInterval(() => {
            progress.value = Tone.Transport.seconds;
            if (Tone.Transport.seconds >= duration) {
                clearInterval(progressTimer);
                isPlaying = false;
                Tone.Transport.stop();
                playPauseBtn.textContent = "▶️ Play DJ";
                progress.value = 0;
            }
        }, 50);
    }
}

function updateDJRoll(notes) {
    const svg = document.getElementById("pianoRoll");
    svg.innerHTML = "";

    notes.slice(0, 100).forEach((note, i) => {
        const rect = document.createElementNS("http://www.w3.org/2000/svg", "rect");
        rect.setAttribute("x", i * 8);
        rect.setAttribute("y", 200 - (note - 40) * 4);
        rect.setAttribute("width", 6);
        rect.setAttribute("height", 4);
        rect.setAttribute("fill", `hsl(${(note % 12) * 30}, 80%, 60%)`);
        svg.appendChild(rect);
    });
}
