#pragma link C++ class vector < float> + ;
#pragma link C++ class vector < vector < float>> + ;
#pragma link C++ class vector < int> + ;
#pragma link C++ class vector < vector < int>> + ;

#ifndef StJetTreeStruct_h
#define StJetTreeStruct_h

const int ConstMax = 5000;
struct StJetTreeStruct
{
    /* There's only one jet per event, due to the structure of our analysis. If an event has more than one jet of interest, that eventid will have a frequency > 1.
    Nothing else should change.
    The event level info we store are            : run id, event id, event refmult, event centrality, event triggers.
    The jet level info we store are              : jetpt, jetcorrectedpt, jeteta, jetphi, jetarea, jetradius, jetenergy, jetnef (neutral energy fraction), fRhoVal used,
                                                   highest energy track pt, number of constituents
    The constituent level info we store are      : trackid, trackpt, tracketa, trackphi, trackpx, trackpy, trackpz
    */

    // Jet Data Info
    int runid;
    int eventid;
    float refmult;
    float grefmult;
    float centrality;
    float refcorr2;
    float mcrefmult;
    float weight;
    vector<unsigned int> triggers;
    vector<double> primaryvertex;
    vector<double> primaryvertexerror;

    float maxtrackpt;
    float maxtowerEt;

    float jetpt;
    float jetcorrectedpt;
    float jeteta;
    float jetphi;
    float jetarea;
    float jetradius;
    float jetenergy;
    float jetnef;
    float fRhoValforjet;
    float fSigValforjet;
    float jethighesttrackpt;
    int numberofconstituents;

    float d0z;
    float d0mass;
    float d0pt;
    float d0phi;
    float d0eta;

    float pionpt;
    float pioneta;
    float pionphi;
    float pioncharge;

    float kaonpt;
    float kaoneta;
    float kaonphi;
    float kaoncharge;

    float lambda_1_half;
    float lambda_1_1;
    float lambda_1_1half;
    float lambda_1_2;
    float lambda_1_3;
    float dispersion;

    float mTrackID[ConstMax];
    float mTrackPt[ConstMax];
    float mTrackEta[ConstMax];
    float mTrackPhi[ConstMax];
    float mTrackPx[ConstMax];
    float mTrackPy[ConstMax];
    float mTrackPz[ConstMax];
    float mTrackCharge[ConstMax];

    void Clear()
    {
        runid = 0;
        eventid = 0;
        refmult = 0;
        grefmult = 0;
        centrality = 0;
        refcorr2 = 0;
        mcrefmult = 0;
        weight = 0;
        triggers.clear();
        primaryvertex.clear();
        primaryvertexerror.clear();

        maxtrackpt = 0;
        maxtowerEt = 0;

        jetpt = 0;
        jetcorrectedpt = 0;
        jeteta = 0;
        jetphi = 0;

        jetarea = 0;
        jetradius = 0;
        jetenergy = 0;
        jetnef = 0;
        fRhoValforjet = 0;
        fSigValforjet = 0;
        jethighesttrackpt = 0;

        numberofconstituents = 0;

        d0z = 0;
        d0mass = 0;
        d0pt = 0;
        d0phi = 0;
        d0eta = 0;

        pionpt = 0;
        pioneta = 0;
        pionphi = 0;
        pioncharge = 0;

        kaonpt = 0;
        kaoneta = 0;
        kaonphi = 0;
        kaoncharge = 0;

        lambda_1_1 = 0;
        lambda_1_1half = 0;
        lambda_1_2 = 0;
        lambda_1_3 = 0;
        dispersion = 0;
    }
};

#endif