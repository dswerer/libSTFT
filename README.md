# dSTFT – A STFT Implementation for Audio Spectrum Analysis

This library provides a convenient C++ interface for computing Short‑Time Fourier Transforms (STFT) of audio files, with flexible data management (disk, memory, or on‑the‑fly) and a uniform `Spectrum` interface. It leverages **libsndfile** for file I/O and **FFTW** for fast Fourier transforms.

---

## Features

- **Multi‑channel support** – process each channel independently.
- **Three load modes**:
  - `FILE` – write spectra and cumulative sums to disk (low memory footprint).
  - `MEMORY` – keep all spectra in RAM for fastest random access.
  - `RUNTIME` – compute FFT on demand (no storage overhead).
- **Uniform `Spectrum` interface** – access frames via `operator[]`, `readNextFrame()`, `seek()`, etc.
- **Cumulative sum (presum)** – each channel provides a prefix‑sum spectrum for fast range‑sum queries.
- **Windowing** – Hann window (configurable, default is rectangular).
- **Optional log‑scale magnitude** – `doLog` flag (default `true`) applies `log2(magnitude + 1)`.
- **High‑frequency boost** – `boostHF` flag linearly amplifies higher bins.
- **Raw PCM access** – retrieve original audio samples through the same `Spectrum` interface.
- **Thread‑safe design** – disk‑based spectra are *not* thread‑safe; use memory or runtime modes for concurrent access, or open separate file handles.

---

## Quick Start

Include the header `stft.h` and link against `libsndfile`, `fftw3`, and `nlohmann/json`.

### 1. Open a sound file
```cpp
#include "stft.h"
using namespace dSTFT;

SoundFile sf("audio.wav");
if (!sf.available) {
    // handle error
}
sf.printInfo();   // display metadata
```

### 2. Create a loader
```cpp
Loader loader(sf, 1024, 2048);   // stride = 1024, frame size = 2048
// Default constructor uses stride=1024, frameSize=2048
// loader.setWindow(hann);       // apply Hann window (optional)
// loader.cdt.doLog = true;      // log‑scale magnitude (default)
// loader.cdt.boostHF = false;   // no HF boost (default)
```

### 3. Generate spectrum files (disk mode)
```cpp
loader.LoadToFile();   // writes: audio.wav-meta.json, audio.wav-c0, audio.wav-c0-ps, ...
```

### 4. Load from disk later
```cpp
SpctCollected sc = Loader::CollectFromFile("audio.wav");
```

### 5. Or load directly into memory
```cpp
SpctCollected sc = loader.LoadToMemory();
```

### 6. Or use runtime (on‑the‑fly) mode
```cpp
SpctCollected sc = loader.LoadRunTime();
// Each access computes FFT from the original file
```

### 7. Access the spectra
```cpp
if (sc.valid) {
    auto& ch0 = sc.channel[0];               // first channel spectrum
    const double* frame0 = (*ch0)[0];        // pointer to first frame
    double bin10 = frame0[10];               // 10th frequency bin

    // Iterate sequentially
    for (uint32_t i = 0; i < sc.meta.frames; ++i) {
        const double* f = ch0->readNextFrame();
        // process...
    }

    // Cumulative sum (presum) – prefix sums over frequency
    auto& ps0 = sc.presum[0];
    double sum_bins_0_to_10 = (*ps0)[0][10]; // sum of first 10 bins of frame 0
}
```

### 8. Retrieve raw PCM samples (original audio)
```cpp
auto raw = loader.getRawPCM(0);   // channel 0
const double* sample = (*raw)[100]; // 100th sample (frame index = sample number)
```

---

## Detailed Class Reference

### `SoundFile`
- Constructors: `SoundFile(const char*)`, copy, move.
- `ReadChunk(start, width)` – reads `width` samples (all channels) from `start`.
- `SplitChannels(raw)` – returns a vector of per‑channel data.
- `getChannel(raw, channelId)` – extracts a single channel.
- `printInfo()` – prints `SF_INFO` fields.

### `Conductor`
- Holds FFTW plan and manages windowing/FFT.
- `doFFT(inputSeries)` – applies window, executes FFT, returns magnitude vector (log or linear).
- `nextFrameFFT()` – reads next chunk, splits channels, calls `doFFT` for each, returns vector of channel spectra.
- `loadMeta()` – returns JSON metadata.
- Members: `doLog` (bool), `boostHF` (bool), `window` (double*), `frameSize`, `stride`, `frameMax`.

### `Loader`
- Constructors: `Loader(SoundFile sf, uint32_t stride=1024, uint32_t frameSize=2048)`.
- `LoadToFile()` – saves spectra and presum to disk.
- `LoadToMemory()` – loads everything into RAM (returns `SpctCollected` with `SpctFromMemory`).
- `LoadRunTime()` – returns `SpctCollected` with `SpctRunTime` (on‑demand FFT).
- `CollectFromFile(fileName)` – static method to read previously saved spectra from disk.
- `Load(LoadMode mode)` – unified entry: `FILE`, `MEMORY`, or `RUNTIME`.
- `setWindow(WindowType)` – currently supports `hann` (others can be added).
- `getRawPCM(channelId)` – returns a `Spectrum*` that provides raw sample values (no FFT).
- `hanning(target, n)` – static helper to fill a Hann window.

### `Spectrum` (abstract base)
- `virtual const double* operator[](int i) = 0` – random access to frame `i`.
- `virtual const double* readNextFrame() = 0` – sequential access.
- `virtual bool seek(uint32_t n) = 0` – set current frame position.
- `virtual uint32_t presentFrame() = 0` – current frame index.
- `virtual metadata meta() = 0` – returns `{binCount, samplerate, frames, stride}`.
- `virtual Spectrum* copyHandle() = 0` – create an independent copy (useful for multi‑threading).

### `SpctCollected`
- Holds vectors of `unique_ptr<Spectrum>` for `channel` and `presum` (cumulative sums).
- `meta` – spectrum metadata.
- `valid` – indicates successful loading.

---

## Important Notes

- **File mode (`FILE`)**:  
  - Generates one binary file per channel for the spectrum (`-cN`) and one for the presum (`-cN-ps`), plus a JSON metadata file.  
  - `SpctFromFile` shares a single `FILE*` internally; **not thread‑safe**. For concurrent reads, use `copyHandle()` to obtain independent handles.

- **Memory mode (`MEMORY`)**:  
  - All spectra are stored as `double*` arrays in `std::vector`.  
  - Fast random access, but memory usage scales with `frames * binCount * channels`.

- **Runtime mode (`RUNTIME`)**:  
  - Computes FFT on every access (`operator[]` or `readNextFrame`).  
  - Low memory, but slower – suitable for interactive visualisation with caching.

- **Presum (cumulative sum)**:  
  - For each frame, `presum[i][j] = sum_{k=0..j} spectrum[i][k]`.  
  - Enables O(1) range‑sum queries (e.g., sum of bins 5–15 = `presum[15] - presum[4]`).

- **Windowing**:  
  - Default is a rectangular window (all ones).  
  - Call `loader.setWindow(hann)` before loading to apply a Hann window.

- **Log‑scale & HF boost**:  
  - `doLog = true` (default) gives `log2(magnitude + 1)` for better visual dynamics.  
  - `boostHF = true` multiplies each bin by `(bin_index / frameSize + 1)` to emphasise high frequencies.

- **Raw PCM**:  
  - `getRawPCM()` returns a `Spectrum` where each “frame” is a single sample (binCount = frameSize).  
  - Useful for waveform rendering or further custom processing.

---

## Example Workflow (All Modes)

```cpp
SoundFile sf("song.wav");
Loader loader(sf, 2048, 4096);
loader.setWindow(hann);
loader.cdt.boostHF = true;   // optional

// Generate disk files
loader.LoadToFile();

// Load from disk (disk mode)
SpctCollected diskSc = Loader::CollectFromFile("song.wav");

// Or load into memory (memory mode)
SpctCollected memSc = loader.LoadToMemory();

// Or use runtime mode
SpctCollected rtSc = loader.LoadRunTime();

// Access first channel, frame 10, bin 5
double val = (*memSc.channel[0])[10][5];
```

---

## Dependencies

- [libsndfile](http://www.mega-nerd.com/libsndfile/) – audio file reading.
- [FFTW](http://www.fftw.org/) – fast Fourier transforms.
- [nlohmann/json](https://github.com/nlohmann/json) – JSON for metadata.

---

## License

This library is part of a sound visualizer project. Feel free to use and adapt it for your own needs.
