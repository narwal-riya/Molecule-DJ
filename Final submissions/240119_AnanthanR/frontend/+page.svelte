<script>
  import * as Tone from 'tone';

  let smiles = '';
  let notes = [];
  let error = '';
  let loading = false;

  async function submitPrompt() {
    error = '';
    notes = [];
    loading = true;

    try {
      const response = await fetch('http://localhost:8000/generate', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ smiles })
      });

      if (!response.ok) {
        error = "Failed to generate music.";
        loading = false;
        return;
      }

      const data = await response.json();
      notes = data.notes;
    } catch (err) {
      error = "Error connecting to backend.";
    } finally {
      loading = false;
    }
  }

  async function checkHealth() {
    try {
      const res = await fetch('http://localhost:8000/health-check');
      const data = await res.json();
      alert(data.status);
    } catch {
      alert("Backend unreachable.");
    }
  }

  function playNotes() {
    const synth = new Tone.Synth().toDestination();
    const now = Tone.now();

    notes.forEach((midi, index) => {
      const note = Tone.Frequency(midi, "midi").toNote();
      synth.triggerAttackRelease(note, "8n", now + index * 0.3);
    });
  }
</script>

<style>
  textarea {
    width: 90%;
    padding: 1rem;
    font-size: 1rem;
  }

  button {
    margin: 0.5rem;
    padding: 0.5rem 1rem;
    font-size: 1rem;
    cursor: pointer;
  }

  .error {
    color: red;
    font-weight: bold;
  }
</style>

<h1>🎵 Molecule DJ</h1>

<textarea
  bind:value={smiles}
  rows="5"
  placeholder="Enter a SMILES string for your molecule...">
</textarea>

<br />

<button on:click={checkHealth}>Check Backend Status</button>
<button on:click={submitPrompt} disabled={loading}>
  {loading ? "Generating..." : "Generate Music"}
</button>

{#if error}
  <p class="error">{error}</p>
{/if}

{#if notes.length > 0}
  <h2>🎼 Generated Notes:</h2>
  <p>{notes.join(", ")}</p>
  <button on:click={playNotes}>▶️ Play</button>
{/if}
