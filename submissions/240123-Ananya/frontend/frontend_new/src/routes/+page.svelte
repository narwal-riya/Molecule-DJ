<script>
	let text = "";
	let audioUrl = "";

	async function checkStatus() {
		const res = await fetch("http://127.0.0.1:8000/health-check");
		const status = await res.text();
		alert(status); // should alert "healthy"
	}

	async function generateAudio() {
		const res = await fetch("http://127.0.0.1:8000/generate", {
			method: "POST",
			headers: {
				"Content-Type": "application/json"
			},
			body: JSON.stringify({ text })
		});

		if (!res.ok) {
			alert("Backend failed with status: " + res.status);
			return;
		}

		const blob = await res.blob();

		// Double-check blob type
		console.log("Blob type:", blob.type);
		console.log("Blob size:", blob.size);

		// Ensure it's served as audio/wav
		const fixedBlob = new Blob([blob], { type: "audio/wav" });
		audioUrl = URL.createObjectURL(fixedBlob);
	}
</script>

<main>
	<h1>Text to Audio Generator</h1>

	<textarea bind:value={text} rows="5" cols="40" placeholder="Enter your text here..."></textarea>

	<div style="margin-top: 1em;">
		<button on:click={checkStatus}>Status</button>
		<button on:click={generateAudio}>Submit</button>
	</div>

	{#if audioUrl}
		<audio controls src={audioUrl}></audio>
	{/if}
</main>

<style>
	main {
		display: flex;
		flex-direction: column;
		align-items: center;
		margin-top: 2em;
		font-family: sans-serif;
	}

	textarea {
		width: 300px;
		padding: 0.5em;
	}

	button {
		margin: 0.5em;
		padding: 0.5em 1em;
		font-size: 1em;
	}
</style>
