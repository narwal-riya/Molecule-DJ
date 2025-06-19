<script>
  let textPrompt = "";
  let audioSrc = "";

  async function checkHealth() {
    try {
      const res = await fetch("http://localhost:8000/health-check", {
        method: "POST",
      });
      const data = await res.json();
      alert("Backend Status: " + data.status);
    } catch (err) {
      alert("Health check failed: " + err);
    }
  }

  async function submitPrompt() {
    try {
      const res = await fetch("http://localhost:8000/generate", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ prompt: textPrompt })
      });

      if (!res.ok) {
        throw new Error("Failed to fetch audio");
      }

      const blob = await res.blob();
      audioSrc = URL.createObjectURL(blob);
    } catch (err) {
      alert("Error generating audio: " + err.message);
    }
  }
</script>

<main>
  <h2>🎵 Molecule DJ</h2>

  <textarea bind:value={textPrompt} rows="5" cols="50" placeholder="Enter music prompt..."></textarea>
  <br /><br />

  <button on:click={checkHealth}>Check Status</button>
  <button on:click={submitPrompt}>Generate</button>

  <br /><br />
  {#if audioSrc}
    <audio src={audioSrc} controls autoplay></audio>
  {/if}
</main>
