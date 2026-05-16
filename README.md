## A STFT implementation for my sound visualizer
## Features:
1. Create spectrum for each channel from raw sound file with few lines(Thanks to FFTW and libsndfile)
2. Manage data in memory or in disk according to your need
3. Uniform interface for both ways above
4. Presum of each spectrum so you can easily compute the sum of a certain range

## How to use:
Just in case if this library is of any help to anyone, here's a quick guide on using it.

---

First, you create a SoundFile object to read your soundfile(include the header btw). This class is actually a wrapped sndfile handle. The SND_INFO struct is accessible with "info".
```c++
    SoundFile sf("path/to/soundfile.wav");
    sf.printInfo();
```
Check sf.available to see if the file is opened correctly.

---

Second, you create a Loader object.
```c++
    Loader ld(sf,1024,2048);
```
This constructor copies the sf object to a internal Conductor object(which wraps fftw). It recreates the SNDFILE handle so that you can still use sf happily. The Conductor is initialized with stride (1024) and frame size (2048) (which are default value).

---

Third, generate spectrum with ld.
```c++
    ld.LoadToFile();
    SpctCollected sc=ld.LoadToMemory();
```
The first line generates 2 files for each channel, one spectrum and one presum (with a '-ps' on its tail). It also generates a json of all the metadata needed.

The second line returns a struct whose definition lies below.
```c++
typedef struct{
    std::vector<std::unique_ptr<Spectrum>> channel;
    std::vector<std::unique_ptr<Spectrum>> presum;
    Spectrum::metadata meta;
    bool valid=false;
}SpctCollected;
```
This is where you fetch the spectrums.

---

To later read a previously saved spectrum from disk:
```c++
    SpctCollected sc=Loader::CollectFromFile("path/to/soundfile.wav");
```

---

The abstracted class Spectrum defines the following functions
```c++
    virtual const double* operator [](int i)=0;
    virtual const double* readNextFrame()=0;
    virtual bool seek(uint32_t n)=0;
    virtual uint32_t presentFrame()=0;
    virtual metadata meta()=0;
```
operator[](int i) is overloaded for random visit to each frame

readNextFrame() and seek() is for higher performance in streaming visits. presentFrame() returns the present frame index.
```c++
    Spectrum* spct=...;
    (*spct)[0][10]; //visiting the 10th bit of frame 0
    double f=spct->readNextFrame(); //now returns a pointer to the first frame
```

---

Altogether, a typical workflow looks like:
```c++
    SoundFile sf("path/to/soundfile.wav");
    Loader ld(sf);
    ld.LoadToFile();
    SpctCollected sc=Loader::CollectFromFiles("path/to/soundfile.wav");
    std::cout<<"10th bin of the first frame of first channels"<<(*sc.channel[0])[0][10];
```

***

Note: The file‑based spectrum (SpctFromFile) is not thread‑safe due to shared internal file pointers. Use LoadToMemory() or open separate file handles per thread if concurrent access is needed.