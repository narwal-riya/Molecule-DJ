<script>
	let inputText = '';
	let audioUrl = '';

	// Check backend health
	async function checkBackend() {
		try {
			const res = await fetch('http://127.0.0.1:12345/health-check', {
				method: 'POST',
			});
			const text = await res.text();
			alert('Backend status: ' + text);
		} catch (err) {
			alert('Error connecting to backend.');
		}
	}

	// Generate audio from text
	async function generateMusic() {
		try {
			const res = await fetch('http://127.0.0.1:12345/generate', {
				method: 'POST',
				headers: { 'Content-Type': 'application/json' },
				body: JSON.stringify({ text: inputText })
			});
			const blob = await res.blob();
			audioUrl = URL.createObjectURL(blob);
		} catch (err) {
			alert('Failed to generate music.');
		}
	}
</script>

<textarea bind:value={inputText} rows="6" cols="40" placeholder="Enter your prompt here..." />

<br />

<button on:click={checkBackend}>Check Backend</button>
<button on:click={generateMusic}>Generate Music</button>

{#if audioUrl}
	<br /><br />
	<audio controls src={audioUrl}></audio>
{/if}
