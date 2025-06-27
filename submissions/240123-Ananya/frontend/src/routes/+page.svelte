<script>
	let text = '';
	let audioUrl = '';

	async function checkStatus() {
		const res = await fetch('http://127.0.0.1:8000/health-check');
		const data = await res.text();
		alert(data); // Should say "healthy"
	}

	async function submitText() {
		console.log("Submitting text:", text);

		const res = await fetch('http://127.0.0.1:8000/generate', {
			method: 'POST',
			headers: {
				'Content-Type': 'application/json'
			},
			body: JSON.stringify({ text })
		});

		const blob = await res.blob();
		console.log("Blob type:", blob.type);
		console.log("Blob size:", blob.size);

		audioUrl = ''; // reset
		const fixedBlob = new Blob([blob], { type: 'audio/wav' });
		audioUrl = URL.createObjectURL(fixedBlob);

		// Optional: auto download for debugging
		const a = document.createElement('a');
		a.href = audioUrl;
		a.download = 'test.wav';
		a.click();
	}
</script>

<textarea bind:value={text} placeholder="Enter text to generate audio..."></textarea>
<br />
<button on:click={checkStatus}>Status</button>
<button on:click={submitText}>Submit</button>

{#if audioUrl}
	<audio controls src={audioUrl}></audio>
{/if}
<script>
	let text = '';
	let audioUrl = '';

	async function checkStatus() {
		const res = await fetch('http://127.0.0.1:8000/health-check');
		const data = await res.text();
		alert(data); // Should say "healthy"
	}

	async function submitText() {
		console.log("Submitting text:", text);

		const res = await fetch('http://127.0.0.1:8000/generate', {
			method: 'POST',
			headers: {
				'Content-Type': 'application/json'
			},
			body: JSON.stringify({ text })
		});

		const blob = await res.blob();
		console.log("Blob type:", blob.type);
		console.log("Blob size:", blob.size);

		audioUrl = ''; // reset
		const fixedBlob = new Blob([blob], { type: 'audio/wav' });
		audioUrl = URL.createObjectURL(fixedBlob);

		// Optional: auto download for debugging
		const a = document.createElement('a');
		a.href = audioUrl;
		a.download = 'test.wav';
		a.click();
	}
</script>

<textarea bind:value={text} placeholder="Enter text to generate audio..."></textarea>
<br />
<button on:click={checkStatus}>Status</button>
<button on:click={submitText}>Submit</button>

{#if audioUrl}
	<audio controls src={audioUrl}></audio>
{/if}
