#include"stft.h"
using namespace dSTFT;
SoundFile::SoundFile(const char* fileName){
    // SF_INFO info;
    this->available=false;
    this->info=(SF_INFO*)calloc(1,sizeof(SF_INFO));
    SNDFILE* sf = sf_open(fileName,SFM_READ,this->info);
    this->fileName.assign(fileName);
    if(!sf){
        std::cerr<<"SOUNDFILE:FAILED TO OPEN SOUND FILE";
        return;
    }
    // std::cout<<"Samplerate:"<<this->info->samplerate<<std::endl;
    // std::cout<<"Channels:"<<this->info->channels<<std::endl;
    // std::cout<<"Frames:"<<this->info->frames<<std::endl;
    
    this->sndfile=sf;
    this->fileName=std::string(fileName);
    this->available=true;
}

SoundFile::SoundFile(const SoundFile& sf){
    this->available=sf.available;
    this->info=(SF_INFO*)calloc(1,sizeof(SF_INFO));
    SNDFILE* sff=sf_open(sf.fileName.c_str(),SFM_READ,this->info);
    if(!sff){
        std::cerr<<"SOUNDFILE:FAILED TO OPEN SOUND FILE";
        return;
    }
    this->fileName=std::string(sf.fileName);
    this->sndfile=sff;
}

SoundFile::~SoundFile(){
    if(info) free(info);
    if(sndfile) sf_close(sndfile);
}

std::vector<double> SoundFile::ReadChunk(uint32_t start,uint32_t width){
    sf_seek(sndfile,start,SEEK_SET);
    // double* buffer=new double[width*this->info->channels];
    std::vector<double> buffer(width*info->channels);
    uint32_t read= sf_read_double(sndfile,buffer.data(),width*info->channels);
    // buffer.resize(read);
    return buffer;
}

std::vector<std::vector<double>> SoundFile::SplitChannels(std::vector<double> rawFrames){
    uint32_t frameCount=rawFrames.size()/info->channels;
    std::vector<std::vector<double>> data(info->channels,std::vector<double>(frameCount,0.));
    for(uint32_t c=0;c<info->channels;c++){
        for(uint32_t i=0;i<frameCount;i++){
            data[c][i]=rawFrames[c+i*info->channels];
        }
    }
    return data;
}

std::vector<double> SoundFile::getChannel(std::vector<double>& rawFrames,uint8_t channelId){
    uint32_t frameCount=rawFrames.size()/info->channels;
    std::vector<double> data(frameCount);
    for(uint32_t i=0;i<frameCount;i++){
        data[i]=rawFrames[channelId+i*info->channels];
    }
    return data;
}

void SoundFile::printInfo(){
    using namespace std;
    cout<<"SoundFile Info-----------------------------\n";
    cout<<"channels - "<<info->channels<<endl;
    cout<<"format - "<<info->format<<endl;
    cout<<"frames - "<<info->frames<<endl;
    cout<<"samplerate - "<<info->samplerate<<endl;
    cout<<"sections - "<<info->sections<<endl;
    cout<<"seekable - "<<info->seekable<<endl;
    cout<<"-------------------------------------------\n";

}

Conductor::Conductor(Conductor& src):
soundFile(src.soundFile),stride(src.stride),frameCount(src.frameCount),frameMax(src.frameMax),frameSize(src.frameSize){
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
    available=true;
    doLog=src.doLog;

    window=(double*)malloc(sizeof(double)*frameSize);
    for(int i=0;i<frameSize;i++){
        window[i]=src.window[i];
    }
}

std::vector<double> Conductor::doFFT(std::vector<double> inputSeries){
    for(uint32_t i=0;i<this->frameSize;i++){
        in[i][0]=inputSeries[i]*window[i];
        in[i][1]=0.;
    }
    fftw_execute(plan);
    std::vector<double> outputSeries;
    outputSeries.reserve(this->frameSize);
    for(uint32_t i=0;i<this->frameSize;i++){
        if(doLog){
            outputSeries.push_back(log2(sqrt(out[i][0]*out[i][0]+out[i][1]*out[i][1])+1.));
        }
        else
        outputSeries.push_back((out[i][0]*out[i][0]+out[i][1]*out[i][1]));
    }
    return outputSeries;
}

std::vector<std::vector<double>> Conductor::nextFrameFFT(){
    auto rawdata=soundFile.ReadChunk(frameCount*stride,frameSize);
    if(rawdata.size()<frameSize*soundFile.info->channels){
        this->end=1;
        // return std::vector<std::vector<double>>();
        if(rawdata.size()==0){
            std::cerr<<"WARNING:EMPTY SAMPLING\n";
            return std::vector<std::vector<double>>();
        }
    }
    int k=rawdata.size();
    rawdata.resize(frameSize*soundFile.info->channels);
    for(int i=k;i<rawdata.size();i++){
        rawdata[i]=0.;
    }

    auto data=soundFile.SplitChannels(std::move(rawdata));
    std::vector<std::vector<double>> output;
    for(auto& channel : data){
        auto phase = doFFT(channel);
        output.push_back(phase);
    }
    frameCount++;
    return output;
}

json Conductor::loadMeta(){
    json meta;
    meta["fileName"]=soundFile.fileName;
    meta["channels"]=soundFile.info->channels;
    meta["frames"]=this->frameMax;
    meta["samplerate"]=soundFile.info->samplerate;
    meta["stride"]=this->stride;
    meta["frameSize"]=this->frameSize;
    meta["binCount"]=frameSize/2+1;
    return meta;
}

Conductor::~Conductor(){
    if(plan) fftw_destroy_plan(plan);
    if(in) fftw_free(in);
    if(out) fftw_free(out);
}

// Loader::Loader(SoundFile&& sf):cdt(std::move(sf),1024,2048){
//     available=false;
//     if(!cdt.available){
//         std::cerr<<"LOADER:FAILED TO CONSTRUCT CONDUCTOR\n";
//         return;
//     }
//     available=true;
// }
// Loader::Loader(SoundFile&& sf,uint32_t stride,uint32_t frameSize):
// cdt(std::move(sf),stride,frameSize){
//     available=false;
//     if(!cdt.available){
//         std::cerr<<"LOADER:FAILED TO CONSTRUCT CONDUCTOR\n";
//         return;
//     }
//     available=true;
// }
Loader::Loader(SoundFile sf,uint32_t stride,uint32_t frameSize):cdt(std::move(sf),stride,frameSize){
    available=false;
    if(!cdt.available){
        std::cerr<<"LOADER:FAILED TO CONSTRUCT CONDUCTOR\n";
        return;
    }
    available=true;
}
Loader::Loader(SoundFile sf):cdt(std::move(sf),1024,2048){
    available=false;
    if(!cdt.available){
        std::cerr<<"LOADER:FAILED TO CONSTRUCT CONDUCTOR\n";
        return;
    }
    available=true;
}

SpctCollected Loader::CollectFromFile(std::string fileName){
    SpctCollected result;
    result.valid = false;

    std::string metaFile = fileName + "-meta.json";
    std::ifstream metaStream(metaFile);
    if (!metaStream.is_open()) {
        std::cerr << "FAILED TO OPEN METAFILE :: " << metaFile << std::endl;
        return result;
    }
    json metaJson;
    metaStream>>metaJson;
    result.meta.binCount = metaJson.at("binCount").get<uint32_t>();
    result.meta.frames   = metaJson.at("frames").get<uint32_t>();
    result.meta.samplerate = metaJson.at("samplerate").get<uint32_t>();
    result.meta.stride=metaJson.at("stride").get<uint32_t>();
    int channels = metaJson.at("channels").get<int>();

    result.channel.reserve(channels);
    result.presum.reserve(channels);

    for (int i = 0; i < channels; ++i) {
        std::string specFile = fileName + "-c" + std::to_string(i);
        std::FILE* fpSpec = fopen(specFile.c_str(), "rb");
        std::string psFile = fileName + "-c" + std::to_string(i) + "-ps";
        std::FILE* fpPs = fopen(psFile.c_str(), "rb");
        if(!fpSpec){
            std::cerr<<"SPECTRUM NOT COMPLETE\n";
            return result;
        }

        result.channel.push_back(std::unique_ptr<SpctFromFile>(new SpctFromFile(fpSpec,result.meta,specFile)));
        if(fpPs) result.presum.push_back(std::make_unique<SpctFromFile>(fpPs, result.meta,psFile));
    }

    result.valid = true;
    return result;
}

int Loader::LoadToFile(){
    json meta=cdt.loadMeta();
    std::ofstream metafile(cdt.soundFile.fileName+"-meta.json");
    metafile<<meta<<std::endl;

    for(int i=0;i<102;i++){
        std::cout<<'-';
    }
    printf("|");
    putchar('\r');
    printf("|");
    
    std::vector<std::FILE*> spectrum;
    std::vector<std::FILE*> presum;
    for(int i=0;i<meta.at("channels");i++){
        spectrum.push_back(fopen((cdt.soundFile.fileName+"-c"+std::to_string(i)).c_str(),"wb"));
        presum.push_back(fopen((cdt.soundFile.fileName+"-c"+std::to_string(i)+"-ps").c_str(),"wb"));        
    }
    uint64_t fc=0;
    uint64_t prog=0;
    while(cdt.end!=1){
        auto data=cdt.nextFrameFFT();
        if(data.size()==0) break;
        
        for(int i=0;i<meta.at("channels");i++){
            fwrite(data[i].data(),sizeof(double),data[i].size()/2+1,spectrum[i]);
            std::vector<double> ps;
            double sum=0.f;
            for(int j=0;j<data[i].size()/2+1;j++){
                sum+=data[i][j];
                ps.push_back(sum);
            }
            fwrite(ps.data(),sizeof(double),ps.size(),presum[i]);
        }
        // printf("%d\n",fc++);
        fc++;
        if(double(fc)/double(int(meta.at("frames")))>double(prog)*0.01){
            putchar('#');
            prog++;
        }

        if(fc>int(meta.at("frames"))){
            // std::cout<<"final frame\n";
            break;
        }
    }

    for(std::FILE* f : spectrum){
        fclose(f);
    }

    putchar('\n');
    return 0;
}

SpctCollected Loader::LoadToMemory(){
    json m__=cdt.loadMeta();
    Spectrum::metadata meta; 
    meta.binCount=m__.at("binCount");
    meta.frames=m__.at("frames");
    meta.samplerate=m__.at("samplerate");
    meta.stride=m__.at("stride");

    SpctCollected o;
    o.meta=meta;

    std::vector<std::vector<double*>>* basec=new std::vector<std::vector<double*>>();
    std::vector<std::vector<double*>>* basep=new std::vector<std::vector<double*>>();

    std::vector<std::vector<double*>>& channels=*basec;
    std::vector<std::vector<double*>>& presum=*basep;
    for(int i=0;i<m__.at("channels");i++){
        channels.push_back(std::vector<double*>());
        presum.push_back(std::vector<double*>());
    }



    for(int i=0;i<102;i++){
        std::cout<<'-';
    }
    printf("|");
    putchar('\r');
    printf("|");
    uint64_t fc=0;
    uint64_t prog=0;

    while(cdt.end!=1){
        auto data=cdt.nextFrameFFT();
        if(data.size()==0) break;

        for(int i=0;i<m__.at("channels");i++){
            double* ch=new double[meta.binCount]();
            double* ps=new double[meta.binCount]();
            memcpy(ch,data[i].data(),meta.binCount*sizeof(double));
            double sum=0.0f;
            for(int j=0;j<meta.binCount;j++){
                sum+=data[i][j];
                ps[j]=sum;
            }
            channels[i].push_back(ch);
            presum[i].push_back(ps);
        }

        fc++;
        if(double(fc)/double(meta.frames)>double(prog)*0.01){
            putchar('#');
            prog++;
        }
        if(fc>meta.frames) break;
    }

    for(int i=0;i<m__.at("channels");i++){
        // SpctFromMemory chi(std::move(channels[i]));
        // SpctFromMemory psi(std::move(presum[i]));
        o.channel.push_back(std::unique_ptr<SpctFromMemory>(new SpctFromMemory(&channels[i],meta)));
        o.presum.push_back(std::unique_ptr<SpctFromMemory>(new SpctFromMemory(&presum[i],meta)));
    }
    o.valid=true;
    putchar('\n');
    return o;
}

SpctCollected Loader::LoadRunTime(){ 
    // SpctCollected o;`
    json m__=cdt.loadMeta();
    Spectrum::metadata meta; 
    meta.binCount=m__.at("binCount");
    meta.frames=m__.at("frames");
    meta.samplerate=m__.at("samplerate");
    meta.stride=m__.at("stride");

    SpctCollected o;
    o.meta=meta;

    for(int i=0;i<cdt.soundFile.info->channels;i++){
        o.channel.push_back(std::unique_ptr<SpctRunTime>(new SpctRunTime(cdt,meta,i)));
        // o.channel[o.channel.size()-1].
        // o.presum.push_back(std::unique_ptr<SpctRunTime>(new SpctRunTime(cdt,meta,i)));
    }
    return o;
}

SpctCollected Loader::Load(LoadMode mode){
    if(mode==LoadMode::FILE){
        LoadToFile();
        return CollectFromFile(cdt.soundFile.fileName);
    }else if(mode==LoadMode::MEMORY){
        return LoadToMemory();
    }else if(mode==LoadMode::RUNTIME){
        return LoadRunTime();
    }

}

SpctFromFile::SpctFromFile(std::FILE* _file,metadata _data,std::string _name){
    this->__meta=_data;
    fp=_file;
    this->name.assign(_name);
    obuf=new double[_data.binCount];
}
SpctFromFile::SpctFromFile(std::FILE* _file,json _data,std::string _name){
    this->__meta.binCount=_data.at("binCount");
    this->__meta.frames=_data.at("frames");
    this->__meta.samplerate=_data.at("samplerate");
    this->__meta.stride=_data.at("stride");
    
    this->name.assign(_name);
    fp=_file;
    obuf=new double[__meta.binCount];
}

Spectrum* SpctFromFile::copyHandle(){
    std::FILE* nfp=fopen(name.c_str(),"rb");
    SpctFromFile* n=new SpctFromFile(nfp,this->__meta,name);
    return n;
}

bool SpctFromFile::seek(uint32_t n){
    frame=n;
    return !fseek(fp,n*__meta.binCount*sizeof(double),SEEK_SET);
}

const double* SpctFromFile::operator[](int i){
    if(!seek(i)){
        readCount=0;
        return nullptr;
    }
    readCount=fread(obuf,sizeof(double),__meta.binCount,fp);
    return this->obuf;
}

const double* SpctFromFile::readNextFrame(){
    readCount=fread(obuf,sizeof(double),__meta.binCount,fp);
    frame++;
    return this->obuf;
}

SpctFromFile::~SpctFromFile(){
    fclose(fp);
    delete[] obuf;
}

uint32_t SpctFromFile::presentFrame(){
    return frame;
}

Spectrum::metadata SpctFromFile::meta(){
    return __meta;
}

SpctFromMemory::SpctFromMemory(std::vector<double*>* _data,metadata meta):
data(_data),frame(0),__meta(meta){}

const double* SpctFromMemory::operator[](int i){
    frame=i;
    return (*data)[i];
}

const double* SpctFromMemory::readNextFrame(){
    return (*data)[frame++];
}

bool SpctFromMemory::seek(uint32_t n){
    frame=n;
    if(frame>(*data).size()) return false;
    return true;
}

uint32_t SpctFromMemory::presentFrame(){
    return frame;
}

SpctFromMemory::~SpctFromMemory(){
    for(auto& d:*data){
        delete[] d;
    }
}

Spectrum* SpctFromMemory::copyHandle(){
    return new SpctFromMemory(data,__meta);
}

Spectrum::metadata SpctFromMemory::meta(){
    return __meta;
}

#define m_PI 3.1415926
void Loader::hanning(double* target,int n){
    if(n<=0) return;
    if(n==1) target[0]=1.;
    double N=double(n-1);
    for(int i=0;i<n;i++){
        target[i]=0.5*(1.0-cos(2.0*m_PI*i/N));
    }
}

void Loader::setWindow(WindowType type){
    if(type==none) return;
    if(type==hann){
        hanning(cdt.window,cdt.frameSize);
    }
}

const double* SpctRunTime::operator[](int i){
    if(i==frame) return obuf;
    auto rawdata=cdt.soundFile.ReadChunk(i*cdt.stride,cdt.frameSize);
    // if(rawdata.size()<cdt.frameSize*cdt.soundFile.info->channels){
        
    // }
    int k=rawdata.size();
    rawdata.resize(cdt.frameSize*cdt.soundFile.info->channels);
    for(int i=k;i<rawdata.size();i++){
        rawdata[i]=0.;
    }

    // auto data=cdt.soundFile.SplitChannels(std::move(rawdata));

    auto ch=cdt.soundFile.getChannel(rawdata,channelId);

    if(mode==ORIGIN){
        memcpy(obuf,ch.data(),ch.size()*sizeof(double));
        return obuf;
    }

    auto fftdata=cdt.doFFT(ch);

    memcpy(obuf,fftdata.data(),__meta.binCount*sizeof(double));
    frame=i;
    return obuf;
}

const double* SpctRunTime::readNextFrame(){ 
    return operator[](++frame);
}

bool SpctRunTime::seek(uint32_t n){
    if(n>cdt.frameCount) return false;
    frame=n;
    return true;
}

uint32_t SpctRunTime::presentFrame(){
    return frame;
}

Spectrum* SpctRunTime::copyHandle(){
    return new SpctRunTime(cdt,__meta,channelId);
}

Spectrum::metadata SpctRunTime::meta(){
    return __meta;
}

SpctRunTime::~SpctRunTime(){
    delete[] obuf;
}

std::unique_ptr<Spectrum> Loader::getRawPCM(uint8_t channelId){
    Spectrum::metadata meta; 
    meta.binCount=cdt.frameSize;
    meta.frames=cdt.frameMax;
    meta.samplerate=cdt.soundFile.info->samplerate;
    meta.stride=cdt.stride;
    SpctRunTime* spct=new SpctRunTime(cdt,meta,channelId);
    spct->mode=SpctRunTime::ORIGIN;
    return std::unique_ptr<Spectrum>(spct);
}