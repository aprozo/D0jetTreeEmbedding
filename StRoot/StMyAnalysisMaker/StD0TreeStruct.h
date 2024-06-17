#pragma link C++ class vector<float> +;
#pragma link C++ class vector<vector<float> >+;
#pragma link C++ class vector<int> +;
#pragma link C++ class vector<vector<int> >+;

#ifndef StD0TreeStruct_h
#define StD0TreeStruct_h

struct StD0TreeStruct
{
    int runid;
    int eventid;
    float refmult;
    float centrality;
    vector <unsigned int> triggers;
    vector <double> primaryvertex;
    vector <double> primaryvertexerror;
    float vzmvpd;

    float d0mass;
    float d0pt;
    float d0eta;
    float d0phi;

    float pionpt;
    float pioneta;
    float pionphi;
    float pioncharge;
    float piondca;
    float nsigmapion;

    float kaonpt;
    float kaoneta;
    float kaonphi;
    float kaoncharge;
    float kaondca;
    float nsigmakaon;
};

#endif