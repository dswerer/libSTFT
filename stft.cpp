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

std::vector<double> Conductor::doFFT(std::vector<double> inputSeries){
    for(uint32_t i=0;i<this->frameSize;i++){
        in[i][0]=inputSeries[i];
        in[i][1]=0.;
    }
    fftw_execute(plan);
    std::vector<double> outputSeries;
    for(uint32_t i=0;i<this->frameSize;i++){
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
    meta["bitCount"]=frameSize/2+1;
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
    result.meta.bitCount = metaJson.at("bitCount").get<uint32_t>();
    result.meta.frames   = metaJson.at("frames").get<uint32_t>();
    result.meta.samplerate = metaJson.at("samplerate").get<uint32_t>();
    int channels = metaJson.at("channels").get<int>();

    result.channel.reserve(channels);
    result.presum.reserve(channels);

    for (int i = 0; i < channels; ++i) {
        std::string specFile = fileName + "-c" + std::to_string(i);
        std::FILE* fpSpec = fopen(specFile.c_str(), "rb");
        std::string psFile = fileName + "-c" + std::to_string(i) + "-ps";
        std::FILE* fpPs = fopen(psFile.c_str(), "rb");
        if(!fpSpec||!fpPs){
            std::cerr<<"SPECTRUM NOT COMPLETE\n";
            return result;
        }

        result.channel.push_back(std::unique_ptr<SpctFromFile>(new SpctFromFile(fpSpec,result.meta)));
        result.presum.push_back(std::make_unique<SpctFromFile>(fpPs, result.meta));
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
    meta.bitCount=m__.at("bitCount");
    meta.frames=m__.at("frames");
    meta.samplerate=m__.at("samplerate");

    SpctCollected o;
    o.meta=meta;

    std::vector<std::vector<double*>> channels;
    std::vector<std::vector<double*>> presum;
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
            double* ch=new double[meta.bitCount]();
            double* ps=new double[meta.bitCount]();
            memcpy(ch,data[i].data(),meta.bitCount*sizeof(double));
            double sum=0.0f;
            for(int j=0;j<meta.bitCount;j++){
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
        o.channel.push_back(std::unique_ptr<SpctFromMemory>(new SpctFromMemory(std::move(channels[i]))));
        o.presum.push_back(std::unique_ptr<SpctFromMemory>(new SpctFromMemory(std::move(presum[i]))));
    }
    o.valid=true;
    putchar('\n');
    return o;
}

SpctCollected Loader::Load(LoadMode mode){
    if(mode==LoadMode::FILE){
        LoadToFile();
        return CollectFromFile(cdt.soundFile.fileName);
    }else if(mode==LoadMode::MEMORY){
        return LoadToMemory();
    }
}

SpctFromFile::SpctFromFile(std::FILE* _file,metadata _data){
    this->__meta=_data;
    fp=_file;
    obuf=new double[_data.bitCount];
}
SpctFromFile::SpctFromFile(std::FILE* _file,json _data){
    this->__meta.bitCount=_data.at("bitCount");
    this->__meta.frames=_data.at("frames");
    this->__meta.samplerate=_data.at("samplerate");
    fp=_file;
    obuf=new double[__meta.bitCount];
}

bool SpctFromFile::seek(uint32_t n){
    frame=n;
    return !fseek(fp,n*__meta.bitCount*sizeof(double),SEEK_SET);
}

const double* SpctFromFile::operator[](int i){
    int n=frame;
    if(!seek(i)){
        readCount=0;
        return nullptr;
    }
    readCount=fread(obuf,sizeof(double),__meta.bitCount,fp);
    seek(n);
    return this->obuf;
}

const double* SpctFromFile::readNextFrame(){
    readCount=fread(obuf,sizeof(double),__meta.bitCount,fp);
    frame++;
    return this->obuf;
}

SpctFromFile::~SpctFromFile(){
    fclose(fp);
    delete obuf;
}

uint32_t SpctFromFile::presentFrame(){
    return frame;
}

Spectrum::metadata SpctFromFile::meta(){
    return __meta;
}

SpctFromMemory::SpctFromMemory(std::vector<double*>&& _data):
data(std::move(_data)),frame(0){}

const double* SpctFromMemory::operator[](int i){
    return data[i];
}

const double* SpctFromMemory::readNextFrame(){
    return data[frame++];
}

bool SpctFromMemory::seek(uint32_t n){
    frame=n;
    if(frame>data.size()) return false;
    return true;
}

uint32_t SpctFromMemory::presentFrame(){
    return frame;
}

SpctFromMemory::~SpctFromMemory(){
    for(auto& d:data){
        delete[] d;
    }
}

Spectrum::metadata SpctFromMemory::meta(){
    return __meta;
}
