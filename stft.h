#pragma once

#include<iostream>
#include<cstdio>
#include<fstream>
#include"Dependencies/include/sndfile.h"
#include"Dependencies/include/fftw3.h"
#include<string>
#include<vector>
#include<cstdint>
#include<utility>
#include<cstring>
#include<memory>
#include<cmath>
#include"Dependencies/include/json.hpp"
using json=nlohmann::json;

namespace dSTFT{
    enum LoadMode{
        FILE,MEMORY
    };
    enum WindowType{
        none,hann
    };
    
    class SoundFile{
        private:
        SNDFILE* sndfile;
        public:
        std::string fileName;
        SF_INFO* info;
        bool available;

        SoundFile(const char*);
        SoundFile(SoundFile&& sfr):fileName(sfr.fileName),info(sfr.info),sndfile(sfr.sndfile){
            sfr.info=nullptr;
            sfr.sndfile=nullptr;
            available=sfr.available;
        };
        SoundFile(const SoundFile& sfr);
        
        //note that start and width refers to the number of 'frames'(or 'samples')
        std::vector<double> ReadChunk(uint32_t start,uint32_t width);

        //this function splits the output of ReadChunk
        std::vector<std::vector<double>> SplitChannels(std::vector<double> rawFrames);
        void printInfo();
        ~SoundFile();
    };

    class Conductor{
        private:
        uint32_t stride;
        uint32_t frameCount;
        uint32_t frameMax;
        fftw_complex *in;
        fftw_complex *out;
        fftw_plan plan;
        public:
        double* window;
        uint32_t frameSize;
        SoundFile soundFile;
        uint8_t end;
        bool available;
        Conductor(SoundFile&& _soundfile,uint32_t _stride,uint32_t _frameSize)
        :soundFile(std::move(_soundfile)),stride(_stride),frameSize(_frameSize),
        frameCount(0),end(0){
            available=false;
            // window=std::vector<double>(frameSize,0.0);
            in=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*frameSize);
            out=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*frameSize);

            if(in==NULL||out==NULL){
                std::cerr<<"CONDUCTOR:INSUFFICIENT AVAILABLE MEMORY\n";
                return;
            }
            if(!soundFile.available){
                std::cerr<<"CONDUCTOR:SOUNDFILE NOT AVAILABLE\n";
                return;
            }

            plan=fftw_plan_dft_1d(frameSize,in,out,FFTW_FORWARD,FFTW_MEASURE);
            if(soundFile.info->frames<frameSize){
                std::cerr<<"WINDOW SIZE BIGGER THAN TOTAL FRAME(SAMPLE) COUNT\n";
                return;
            }
            frameMax=(soundFile.info->frames-frameSize)/stride;
            available=true;

            window=new double[_frameSize];
            for(int i=0;i<_frameSize;i++){
                window[i]=1.;
            }
        }
        std::vector<double> doFFT(std::vector<double> inputSeries);
        std::vector<std::vector<double>> nextFrameFFT();
        json loadMeta();
        ~Conductor();
    };

    class Spectrum{
        public:
        typedef struct {
            uint32_t binCount;
            uint32_t samplerate;
            uint32_t frames;
            uint32_t stride;
        }metadata;
        virtual const double* operator [](int i)=0;
        virtual const double* readNextFrame()=0;
        virtual bool seek(uint32_t n)=0;
        virtual uint32_t presentFrame()=0;
        virtual metadata meta()=0;
        virtual ~Spectrum(){}
    };

    typedef struct{
        std::vector<std::unique_ptr<Spectrum>> channel;
        std::vector<std::unique_ptr<Spectrum>> presum;
        Spectrum::metadata meta;
        bool valid=false;
    }SpctCollected;

    class SpctFromFile:public Spectrum{
        public:
        const double* operator[](int i) override;
        const double* readNextFrame() override;
        bool seek(uint32_t n) override;
        uint32_t presentFrame() override;
        metadata meta() override;
        SpctFromFile(std::FILE* _file,json _meta);
        SpctFromFile(std::FILE* _file,metadata _meta);
        ~SpctFromFile() override;
        private:
        metadata __meta;
        uint32_t frame;
        size_t readCount;
        std::FILE* fp;
        std::FILE* fprandom;
        double* obuf;
    };

    class SpctFromMemory:public Spectrum{
        public:
        const double* operator[](int i) override;
        const double* readNextFrame() override;
        bool seek(uint32_t n) override;
        uint32_t presentFrame() override;
        metadata meta() override;
        SpctFromMemory(std::vector<double*>&& data,metadata meta);
        ~SpctFromMemory() override;
        private:
        std::vector<double*> data;
        uint32_t frame;
        metadata __meta;
    };


    class Loader{
        private:
        Conductor cdt;
        public:
        bool available;
        // Loader(SoundFile&& sf);
        // Loader(SoundFile&& sf,uint32_t stride,uint32_t frameSize);
        
        Loader(SoundFile sf,uint32_t stride,uint32_t frameSize);
        Loader(SoundFile sf);

        int LoadToFile();
        SpctCollected LoadToMemory();
        static SpctCollected CollectFromFile(std::string fileName);
        SpctCollected Load(LoadMode mode);
        void setWindow(WindowType type);

        static void hanning(double* target,int n);
        // ~Loader();
    };
}