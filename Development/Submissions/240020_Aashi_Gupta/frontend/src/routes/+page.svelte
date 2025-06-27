<script>
    import "../style.css";
    let userInput = "";
    let audioSrc = "";

    async function checkHealth() {
        const res = await fetch("http://127.0.0.5:8000/health-check");
        const text = await res.text();
        alert(text);
    }

    async function submitText() {
        const res = await fetch("http://127.0.0.5:8000/generate", {
            method: "POST",
            headers: {
                "Content-Type": "application/json"
            },
            body: JSON.stringify({ text: userInput })
        });
        const blob = await res.blob();
        audioSrc = URL.createObjectURL(blob);
    }
</script>

<h1>MolDJ UI</h1>
<textarea bind:value={userInput} placeholder="Type something..."></textarea>
<br>
<button on:click={checkHealth}>Check Health</button>
<button on:click={submitText}>Generate</button>
<br>
<audio controls src={audioSrc}></audio>
