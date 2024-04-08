#include "MCSystem.h"
#include <cmath>
#include <string>
#include <fstream>
#define _USE_MATH_DEFINES
using MCSys::MCSystem;

namespace Alloy
{
    

        //----------------------------------------------------------------------------------------------------------

    class Alloy : public MCSystem
    {
    public:

        //----------------------------------------------------------------------------------------------------------

        Eigen::RowVectorXd spin;
        double spinI;
        double dspin;
        double dspinA;
        double sig;
        double J;

        //----------------------------------------------------------------------------------------------------------

        std::vector<long double> E;
        std::vector<long double> C;
        std::vector<std::vector<long double>> g;


        bool getE;
        std::vector<long double> Ee;
        std::vector<long double> Ce;
        std::vector<std::vector<long double>> ge;

        bool getAC;
        std::vector<long double> EAC;
        std::vector<long double> CAC;
        std::vector<std::vector<long double>> gAC;

        std::vector<long double> E1;
        int saveCount;
        
        
        std::vector<std::vector<long double>> _si;
        std::vector<std::vector<long double>> _sj;
        std::vector<std::vector<long double>> _sisj;



        std::vector<long double> gBins;
        std::vector<double> gBinsVol;
        double gBinDel;
        bool gBinsGen;

        //----------------------------------------------------------------------------------------------------------

        Eigen::RowVectorXd m_Interactions;
        bool m_Int;
        std::vector<double> SpinOptions;
        std::vector<double> SpinOptionsDist;
        int maxSpin;
        int minSpin;

        
        int targetIndex;
        Eigen::RowVectorXd targetSq;
        double targetSpin;
        double effectiveRadii;
        //----------------------------------------------------------------------------------------------------------

        Alloy(Eigen::MatrixXd _points, double _J = 1, double _sig = std::numeric_limits<double>::infinity(), bool _PBCImaging = true, Eigen::VectorXd _PBCSpacing = Eigen::VectorXd::Zero(0), std::vector<double> _SpinOptions = {-1,1}, std::vector<double> _SpinOptionsDist = {0.5,0.5}, int splitCount = 2, Eigen::MatrixXd tMatrix = Eigen::RowVectorXd::Zero(0), bool classic = false, float _effectiveRadii = 0., bool stdev=false,bool ac = false)
        {  
            getE = stdev;
            getAC = stdev && ac;
            effectiveRadii = _effectiveRadii;
            SpinOptions = _SpinOptions;
            SpinOptionsDist = _SpinOptionsDist;
            double SpinDistSum = std::accumulate(SpinOptionsDist.begin(), SpinOptionsDist.end(), 0.0);
            for (int i = 0; i < (int) SpinOptionsDist.size(); i++)
                SpinOptionsDist[i] = SpinOptionsDist[i]/SpinDistSum;
            
            

            J = _J;
            maxSpin = 0;
            minSpin = RAND_MAX;
            for (double x : SpinOptions)
            {
                if (maxSpin < std::abs(x))
                    maxSpin = std::abs(x);
                if (minSpin > std::abs(x))
                    minSpin = std::abs(x);
            }

            if (!(tMatrix.cols() == tMatrix.rows()) || !(tMatrix.cols() == _points.rows()))
                tMatrix = Eigen::MatrixXd::Identity(_points.rows(),_points.rows());

            if (!(_PBCSpacing.cols() == 1 && _PBCSpacing.rows() == _points.rows())) {_PBCSpacing = Eigen::VectorXd::Constant(_points.rows(),1);}

            SplitInitially = !classic;

            if (!classic)
            {
            DepthFunctions.push_back([this](int idx) -> Eigen::Vector2d { return this->EnergyBounds0(idx); });
            DepthFunctions.push_back([this](int idx) -> Eigen::Vector2d { return this->EnergyBounds1(idx); });
            DepthFunctions.push_back([this](int idx) -> Eigen::Vector2d { return this->EnergyBounds2(idx); });
            }
            DepthFunctions.push_back([this](int idx) -> Eigen::Vector2d { return this->EnergyBounds3(idx); });


            int _recDepth = DepthFunctions.size();
            if (SplitInitially)
                _recDepth++;


            BuildPoints(_points, _PBCImaging, _PBCSpacing, _recDepth, splitCount, tMatrix);
            sig = _sig;


            E1 = {};
            g = {};
            _si = {};
            _sj = {};
            _sisj = {};
            gBins = {};
            gBinsVol = {};
            gBinsGen = false;
            m_Int = false;
            m_Interactions = Eigen::RowVectorXd::Zero(_points.cols());
        }

        int RandSpin()
        {
            std::vector<int> banned = {};
            int choice = 0;
            double bannedDist = 0;
            for (int i = 0; i < (int) SpinOptions.size(); i++)
            {
                if (!((double) (spin.array() == SpinOptions[i]).count() < SpinOptionsDist[i] * (double) spin.cols()))
                {
                    banned.push_back(i);
                    bannedDist += SpinOptionsDist[choice];
                    if (i == choice) {choice++;}
                }
            }

            double seed = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
            seed -= bannedDist;

            while (choice < (int) SpinOptions.size())
            {
                if (choice < (int) SpinOptions.size())
                    seed -= SpinOptionsDist[choice];
                else
                    break;
                if (seed <= 0 || choice + 1 == (int) SpinOptions.size())
                    return SpinOptions[choice];
                bool isBanned = true;
                while (isBanned)
                {
                    choice++;
                    isBanned = false;
                    for (int x:banned)
                        if (x == choice) isBanned = true;
                }

            }
            return SpinOptions[0];
        }

        virtual void ResetBody(bool systemReset = true) override
        {
            E1 = {};
            _si = {};
            _sj = {};
            _sisj = {};
            saveCount = 0;
            if (systemReset)
            {
                Eigen::RowVectorXi spinIdx = Eigen::RowVectorXi::Zero(points._points.cols());
                {
                Eigen::RowVectorXd _spinIdx = (Eigen::RowVectorXd::Random(points._points.cols()).array()).abs();
                double _count = 0.;
                for (int i = 0; i < (int) SpinOptionsDist.size()-1;i++)
                {
                    _count += SpinOptionsDist[i];
                    spinIdx = (_spinIdx.array() >= _count).select(i+1,spinIdx);
                }
                }
                
                for (int i = 0; i < (int) SpinOptionsDist.size()-1;i++)
                {
                    int n = std::ceil((spinIdx.array() == i).count() - SpinOptionsDist[i] * (double) points._points.cols());
                    if (n > 0)
                    {
                        Eigen::RowVectorXd spinIdxN = (spinIdx.array() == i).select(Eigen::RowVectorXd::Random(points._points.cols()),std::numeric_limits<double>::infinity());
                        while (n > 0)
                        {
                            
                            Eigen::RowVectorXd::Index minIndex;
                            spinIdxN.minCoeff(&minIndex);
                            spinIdx(minIndex) = spinIdx(minIndex)+1;
                            spinIdxN(minIndex)=std::numeric_limits<double>::infinity();
                            n--;
                        }
                    }
                }
                for (int i = (int) SpinOptionsDist.size()-1; i > 0;i--)
                {
                    int n = std::ceil((spinIdx.array() == i).count() - SpinOptionsDist[i] * (double) points._points.cols());
                    if (n > 0)
                    {
                        Eigen::RowVectorXd spinIdxN = (spinIdx.array() == i).select(Eigen::RowVectorXd::Random(points._points.cols()),std::numeric_limits<double>::infinity());
                        while (n > 0)
                        {
                            Eigen::RowVectorXd::Index minIndex;
                            spinIdxN.minCoeff(&minIndex);
                            spinIdx(minIndex) = spinIdx(minIndex)-1;
                            spinIdxN(minIndex)=std::numeric_limits<double>::infinity();
                            n--;
                        }
                    }
                }
                
                spin = Eigen::RowVectorXd::Zero(points._points.cols());
                for (int i = 0; i < (int) points._points.cols(); i++)
                {
                    spin(i) = SpinOptions[spinIdx(i)];
                }
            }
        }

        //----------------------------------------------------------------------------------------------------------

        virtual void OnSuccess() override
        {
            spin(0, pointIdx) = -spin(0, pointIdx);
            spin(0, targetIndex) = -spin(0, targetIndex);
        }

        virtual void OnFail() override
        {
            return;
        }

        //----------------------------------------------------------------------------------------------------------

        void Initialise() override
        {
            spinI = spin(0, pointIdx);
            m_Int = false;
            m_Interactions *= 0;

            Eigen::RowVectorXd _Sq = (spin.array() == spinI).select(std::numeric_limits<double>::infinity(),points.Sq);
            Eigen::RowVectorXd::Index minIndex;
            double _ = _Sq.minCoeff(&minIndex);

            _Sq = (_Sq.array() > _*2).select(std::numeric_limits<double>::infinity(),Eigen::RowVectorXd::Random(_Sq.cols()));
            _Sq.minCoeff(&minIndex);

            targetIndex = minIndex;
            targetSpin = spin(0,targetIndex);
            
            if (targetSpin == spinI)
            {
                int _pointIdx = (int) round((static_cast<double>(rand()) / static_cast<double>(RAND_MAX)) * (points._points.cols()-1));
                pointIdx = _pointIdx;
                points.Initialise(pointIdx);
                Initialise();
                return;
            }

            dspin = 1 - targetSpin/spinI;
            dspinA = std::abs(dspin);

            targetSq = points.SqPBCImageOverride(points._points,targetIndex);
            targetSq(0,targetIndex) = std::numeric_limits<double>::infinity();
            points.Sq(0,targetIndex) = std::numeric_limits<double>::infinity();
            targetSq(0,pointIdx) = std::numeric_limits<double>::infinity();
        }

        void UpdateSplitMatrices(int idx) override
        {
            return;
        }

        //----------------------------------------------------------------------------------------------------------
        
        /*
        Eigen::RowVectorXd mInteraction;
        virtual void PreCalculate() override
        {
            mInteraction = precalcInteraction();
        }

        //----------------------------------------------------------------------------------------------------------

        Eigen::RowVectorXd precalcInteraction()
        {
            double p = -0.5*((double)points._points.rows() + sig);
            Eigen::RowVectorXd interactions = (((points._points).array().square().colwise().sum()).pow(p)).unaryExpr([](double v) { return std::isfinite(v)? v : 0.0; });
            return interactions;
        }

        //----------------------------------------------------------------------------------------------------------
    
        Eigen::Vector2d mEnergyBounds0(int idx)
        {
            Eigen::Vector2d out;
            double _dE = 2 * ((double)mSize(mappedIndex(0,idx))) * MaxContributionInteraction(idx);
            out << -_dE, _dE;
            return out;
        }

        Eigen::Vector2d mEnergyBounds1(int idx)
        {
            Eigen::Vector2d spins = mSpinCount(all,mappedIndex(0,idx)).cast<double>();
            if (!spinI) {spins.reverseInPlace();}
            double minInt = MinContributionInteraction(idx);
            double maxInt = MaxContributionInteraction(idx);
            Eigen::Matrix2d IntMul;
            IntMul << minInt,-maxInt,
                      maxInt,-minInt;
            //double dEmin = 2 * (spins((1 + spinI) / 2, 0) * minInt - spins((1 - spinI) / 2, 0) * maxInt);
            //double dEmax = 2 * (spins((1 + spinI) / 2, 0) * maxInt - spins((1 - spinI) / 2, 0) * minInt);
            return 2*IntMul*spins;
        }

        Eigen::Vector2d mEnergyBounds2(int idx)
        {
            Eigen::Vector2d out;
            Eigen::Vector2d spins = mSpinCount(all,mappedIndex(0,idx)).cast<double>();
            if (!spinI) {spins.reverseInPlace();}
            double totInt = mInteraction(mappedIndex(0,idx));
            double minInt = MinContributionInteraction(idx);
            double maxInt = MaxContributionInteraction(idx);
            double dEmin = 2 * std::max(totInt - 2 * spins(1, 0) * maxInt, 2 * spins(0, 0) * minInt - totInt);
            double dEmax = 2 * std::min(totInt - 2 * spins(1, 0) * minInt, 2 * spins(0, 0) * maxInt - totInt);
            out << dEmin, dEmax;
            return out;
        }
        */    

        Eigen::Vector2d EnergyBounds0(int idx)
        {
            Eigen::Vector2d out;
            Eigen::Vector2d maxminInt = MaxMinContributionInteraction(idx);
            double maxInt = std::max((maxminInt(0)),-(maxminInt(1)));
            double _dE = std::abs(J) *std::abs(spinI)*dspinA* ((double)Size(idx)) * maxInt;
            out << -_dE, _dE;
            
            return out;
        }

        Eigen::Vector2d EnergyBounds1(int idx)
        {
            Eigen::Vector2d spins = SpinCount(idx).cast<double>();
            Eigen::Vector2d maxminInt = MaxMinContributionInteraction(idx);
            double maxInt = maxminInt(0);
            double minInt = maxminInt(1);
            Eigen::Matrix2d IntMul;
            IntMul << minInt,-maxInt,
                      maxInt,-minInt;
            return std::abs(J)*dspinA*std::abs(spinI)*IntMul*spins;
        }

        Eigen::Vector2d EnergyBounds2(int idx)
        {
            Eigen::Vector2d out;
            Eigen::Vector2d spins = SpinCount(idx).cast<double>();
            double totInt = Interaction(idx);
            Eigen::Vector2d maxminInt = MaxMinContributionInteraction(idx);
            double maxInt = maxminInt(0);
            double minInt = maxminInt(1);

            double dEmin = std::abs(spinI)*std::abs(J)*dspinA * std::max(minSpin*totInt - 2 * spins(1, 0) * maxInt
                                        ,-maxSpin*totInt + 2 * spins(0, 0) * minInt);

            double dEmax = std::abs(spinI)*std::abs(J)*dspinA * std::min(maxSpin*totInt - 2 * spins(1, 0) * minInt
                                        ,-minSpin*totInt + 2 * spins(0, 0) * maxInt);
            out << dEmin, dEmax;
            
            return out;
        }

        Eigen::Vector2d EnergyBounds3(int idx)
        {
            return Eigen::Vector2d::Constant(EnergyOfFlip(idx));
        }

        //----------------------------------------------------------------------------------------------------------

        int Size(int idx)
        {
            int size = std::get<1>(points.slice(0,idx));
            int cutoff = std::get<0>(points.slice(0,idx));
            if (pointIdx >= cutoff && pointIdx - size < cutoff)
                size -= 1;
            if (targetIndex >= cutoff && targetIndex - size < cutoff)
                size -= 1;
            return size;
        }

        Eigen::Vector2i SpinCount(int idx)
        {
            Eigen::Vector2i out;
            int size = std::get<1>(points.slice(0,idx));
            int cutoff = std::get<0>(points.slice(0,idx));
            Eigen::RowVectorXd spinIdx = spin.block(0,cutoff,spin.rows(),size);
            int count = (spinIdx.array()*spinI > 0).count();
            out << count, size-count;
            if (pointIdx >= cutoff && pointIdx - size < cutoff)
                out(0,0) -= 1;
            if (targetIndex >= cutoff && targetIndex - size < cutoff)
                out(1,0) -= 1;
            return out;
        }

        Eigen::Vector2d MaxMinContributionInteraction(int idx, bool keepsign = true)
        {
            double p = -0.5*((double)points._points.rows() + sig);
            double interactionMax;
            double interactionMin;
            interactionMax = maxSpin*std::pow(points.getSmallestDistSq(idx),p)-minSpin*std::pow(targetgetLargestDistSq(idx),p);
            interactionMin = minSpin*std::pow(points.getLargestDistSq(idx),p)-maxSpin*std::pow(targetgetSmallestDistSq(idx),p);
            Eigen::Vector2d out;
            out << interactionMax , interactionMin;
            return out;
        }

        double Interaction(int idx)
        {
            double p = -0.5*((double)points._points.rows() + sig);
            Eigen::RowVectorXd interactions = (points.getDistanceSq(idx).array().pow(p) - targetgetDistanceSq(idx).array().pow(p)).unaryExpr([](double v) { return std::isfinite(v)? v : 0.0; });

            m_Int = true;
            m_Interactions.block(0,std::get<0>(points.slice(0,idx)),spin.rows(),std::get<1>(points.slice(0,idx))) = interactions;

            return interactions.sum();
        }

        double EnergyOfFlip(int idx)
        {
            double p = -0.5*((double)points._points.rows() + sig);
            Eigen::RowVectorXd interactions;
            Eigen::RowVectorXd spinArr = spin.block(0,std::get<0>(points.slice(0,idx)),spin.rows(),std::get<1>(points.slice(0,idx)));
            if (!m_Int)
            {
                interactions = (points.getDistanceSq(idx).array().pow(p) - targetgetDistanceSq(idx).array().pow(p)).unaryExpr([](double v) { return std::isfinite(v)? v : 0.0; });
            }
            else
                interactions =  m_Interactions.block(0,std::get<0>(points.slice(0,idx)),spin.rows(),std::get<1>(points.slice(0,idx)));
            double tsinteractions = (spinArr.array() * interactions.array()).sum();
            
            return spinI*J*dspin*tsinteractions;
        }

        //----------------------------------------------------------------------------------------------------------

        virtual void SaveInfo(double T) override
        {
            long double n1 = 1./((long double)(points._points.cols()));

            long double Em = mean(E1);

            std::function<long double (std::vector<long double>,std::vector<long double>)> Cfun = [this,T](std::vector<long double> V1, std::vector<long double> V2) -> long double { return (this->mean(V1) - std::pow(this->mean(V2),2))/(T*T); };
            std::function<long double (std::vector<long double>,std::vector<long double>,std::vector<long double>)> gfun = [this,T](std::vector<long double> V1, std::vector<long double> V2, std::vector<long double> V3) -> long double { return this->mean(V1) - this->mean(V2)*this->mean(V3); };

            std::vector<long double> E2_ = prod(E1,E1);

            long double Cm = Cfun(E2_,E1);

            C.push_back(Cm*n1);
            E.push_back(Em*n1);

            if (getE)
            {
                long double Eerr = StDev(E1,Em);

                std::vector<long double> CJK = JackKnife(Cfun,E2_,E1);
                long double Cerr = StDev(CJK,Cm);

                Ee.push_back(Eerr*n1);

                Ce.push_back(((long double) (CJK.size()-1) / std::sqrt((long double) CJK.size())) *Cerr*n1);

                if (getAC)
                {
                    EAC.push_back(autocorrelation(E1,Em,Eerr,1));

                    CAC.push_back(autocorrelation(CJK,Cm,Cerr,1));
                }
            }

            std::vector<long double>  __g = {};
            std::vector<long double>  __ge = {};
            std::vector<long double>  __gAC = {};
            for (int i = 0; i < (int) _si.size(); i++)
            {
                long double gm = gfun(_sisj[i],_si[i],_sj[i]);
                __g.push_back(gm);
                if (getE)
                {
                    std::vector<long double> gJK = JackKnife3(gfun,_sisj[i],_si[i],_sj[i]);
                    long double gerr = StDev(gJK,gm);
                    __ge.push_back(gerr);
                    if (getAC)
                    {
                        __gAC.push_back(autocorrelation(gJK,gm,gerr,1));
                    }
                }
            } 
            g.push_back(__g);
            if (getE)
            {
                ge.push_back(__ge);
                if (getAC)
                {
                    gAC.push_back(__gAC);
                }
            }
            _si = {};
            _sj = {};
            _sisj = {};
            E1 = {};
            saveCount = 0;
        }
    
        double Volume(double rmax,double rmin)
        {
            if (points._points.rows() % 2 == 0)
            {
                int k = points._points.rows()/2;
                double coeff = std::pow(3.1415926535898,k)/((double) Factorial(k));
                double delta = (std::pow(std::abs(rmax),2*k-1) + std::pow(std::min(rmax,std::abs(rmin)),2*k+1-1))/2;
                return points._points.rows()*coeff*delta;
            }
            else
            {
                int k = (points._points.rows()-1)/2;
                double coeff = 2.*std::pow(4.*3.1415926535898,k)*(double) Factorial(k)/((double) Factorial(2*k+1));
                double delta = (std::pow(std::abs(rmax),2*k-1) + std::pow(std::min(rmax,std::abs(rmin)),2*k+1-1))/2;
                return points._points.rows()*coeff*delta;

            }
        }
        int Factorial(int n)
        {
            int fact = 1;
            for(int i = 1; i <= n; ++i) {
                fact *= i;
            }
            return fact;
        }

        virtual void StoreInfo() override
        {
            std::tuple<long double,std::vector<long double>,std::vector<long double>,std::vector<long double>> Ene_sisj = GetInfo();
            long double Ene = std::get<0>(Ene_sisj);
            std::vector<long double> __si = std::get<1>(Ene_sisj);
            std::vector<long double> __sj = std::get<2>(Ene_sisj);
            std::vector<long double> __sisj = std::get<3>(Ene_sisj);
            E1.push_back(Ene);
            if (_si.size() < gBins.size())
            {
                for (int i = 0; i < (int) gBins.size(); i++)
                {
                    _si.push_back({__si[i]});
                    _sj.push_back({__sj[i]});
                    _sisj.push_back({__sisj[i]});
                }
            }
            else
            {
                for (int i = 0; i < (int) gBins.size(); i++)
                {
                    _si[i].push_back(__si[i]);
                    _sj[i].push_back(__sj[i]);
                    _sisj[i].push_back(__sisj[i]);
                }
            }
            saveCount++;
            return;
        }

        std::tuple<long double,std::vector<long double>,std::vector<long double>,std::vector<long double>> GetInfo ()
        {
            std::vector<long double> __si = {};
            std::vector<long double> __sj = {};
            std::vector<long double> __sisj = {};
            std::vector<int> countN = {};
            long double TotEn = 0;
            std::tuple<long double,std::vector<long double>,std::vector<long double>,std::vector<long double>> out;
            for (int i = 0; i < points._points.cols(); i++){
                double spini = spin(0,i);
                Eigen::MatrixXd Sq = points.SqPBCImageOverride(points._points,i);
                double p = -0.5*((double)points._points.rows() + sig);
                Eigen::RowVectorXd interactions = spin.array() * (Sq.array().pow(p)).unaryExpr([](double v) { return std::isfinite(v)? v : 0.0; });
                double tsinteractions = spini * interactions.sum();
                TotEn = TotEn - (long double) J*tsinteractions;

                Sq = Sq.array().sqrt();
                if (!gBinsGen)
                {
                    gBins = GenerategBins(0.,std::pow((double) points._points.cols(),1/((double)points._points.rows()))/2,4*(int)std::ceil(std::pow((double) points._points.cols(),1/((double)points._points.rows()))));
                    //gBins = GenerategBins(Sq.minCoeff(),Sq.maxCoeff(),4*(int)std::ceil(std::sqrt((double) points._points.rows())*std::pow((double) points._points.cols(),1/((double)points._points.rows()))));
                    /*
                    for (int j = 0; j < (int) gBins.size(); j++){
                        gBinsVol.push_back(Volume(gBins[j]+gBinDel,gBins[j]-gBinDel));
                    }
                    */
                    effectiveRadii = std::max(effectiveRadii,gBinDel);
                    gBinsGen = true;
                }
                if (__si.size() < gBins.size())
                {
                    for (int j = 0; j < (int) gBins.size(); j++){
                        __si.push_back(spini);
                        double _N = ((Sq.array() - gBins[j]).abs() <= effectiveRadii ).count();
                        countN.push_back(0);
                        if (_N == 0)
                        {
                            __sj.push_back(0);
                            __sisj.push_back(0);
                        }
                        else
                        {
                            double sjC = (((Sq.array() - gBins[j]).abs() <= effectiveRadii ).select(spin,0).sum());
                            __sj.push_back(sjC);
                            __sisj.push_back(sjC*spini);
                        }
                    }
                }
                else
                {
                    for (int j = 0; j < (int) gBins.size(); j++){
                        double _N = ((Sq.array() - gBins[j]).abs() <= effectiveRadii ).count();
                        __si[j] = __si[j] + (spini);
                        countN[j] = countN[j] + (_N);
                        if (!(_N == 0))
                        {
                            double sjC = (((Sq.array() - gBins[j]).abs() <= effectiveRadii ).select(spin,0).sum());
                            __sj[j] = __sj[j] + (sjC);
                            __sisj[j] = __sisj[j] + (sjC*spini);
                        }
                    }
                }
            }
            for (int j = 0; j < (int) gBins.size(); j++){
                __si[j] = __si[j]/ ((double)points._points.cols());
                __sj[j] = __sj[j]/ (double) countN[j];
                __sisj[j] = __sisj[j]/ (double) countN[j];
            }
            out = {TotEn,__si,__sj,__sisj};
            return out;
        }

        std::vector<long double> GenerategBins(double _min,double _max,int _size)
        {
            std::vector<long double> linspaced;
            gBinDel = ((_max-_min)/(_size-1))/2;
            long double start = static_cast<long double>(_min);
            long double end = static_cast<long double>(_max+2*gBinDel);
            long double num = static_cast<long double>(_size+1);

            if (num == 0) { return linspaced; }
            if (num == 1) 
                {
                linspaced.push_back(start);
                return linspaced;
                }

            long double delta = (end - start) / (num - 1);

            for(int i=0; i < num-1; ++i)
                {
                linspaced.push_back(start + delta * i);
                }
            linspaced.push_back(end); 
            return linspaced;
        }

        /*
        std::vector<int> RadDistFun()
        {
            std::vector<int> __sisj = {};
            for (int i = 0; i < points._points.cols(); i++){
                Eigen::MatrixXd Sq = points.SqPBCImageOverride(points._points,i);
                if (!gBinsGen)
                {
                    gBins = GenerategBins(Sq.minCoeff(),Sq.maxCoeff(),std::sqrt(points._points.cols()));
                    gBinsGen = true;
                }
                for (int j = 0; j < (int) gBins.size(); j++){
                    __sisj.push_back(((Sq.array() - gBins[j]).abs() <= gBinDel ).count());
                }
            }

        }
        
        double TotalEnergy()
        {
            double TotEn = 0;
            for (int i = 0; i < points._points.cols(); i++){
                Eigen::MatrixXd Sq = points.SqPBCImageOverride(points._points,i);
                double p = -0.5*((double)points._points.rows() + sig);
                Eigen::RowVectorXd interactions = spin.array() * (Sq.array().pow(p)).unaryExpr([](double v) { return std::isfinite(v)? v : 0.0; });
                double tsinteractions = spin(0,i) * interactions.sum();
                TotEn -= tsinteractions;
            }
            return TotEn/2;
        }
        */

        //----------------------------------------------------------------------------------------------------------

        void WriteState(std::string file, double T) override
        {
            std::ofstream outdata(file + "_" + std::to_string(T) + "_" + std::to_string(CycleCount) +".csv");

            if ( outdata.is_open() ) 
            { 
                int i;
                int j;
                outdata << "Ndim , " << points._points.rows() << std::endl;
                outdata << "Temperature , " << T << std::endl;
                outdata << "Cycle Number , " << CycleCount << std::endl;
                outdata << "N , " << points._points.cols() << std::endl;
                for (j=0; j<points._points.rows(); ++j)
                {
                    outdata << "X" << j << "max" << " , " << points.PBCBarrier(j,0) << std::endl;
                }
                for (i=0; i<points._points.rows(); ++i)
                {
                    for (j=0; j<points._points.rows(); ++j)
                    {
                        if (j > 0) {outdata << " , ";}
                        outdata << points.tMatrix(i,j);
                    }
                    outdata << std::endl;
                }

                for (int i = 0; i < (int) SpinOptions.size(); i++)
                {
                    if (i > 0) {outdata << " , ";}
                    outdata << SpinOptions[i];
                }
                outdata << std::endl;

                Eigen::MatrixXd _points = points._points;


                for (j=0; j<points._points.rows(); ++j)
                {
                    if (j > 0) {outdata << " , ";}
                    outdata << "X" << j+1;
                }
                outdata << " , " << "Spin";
                outdata << std::endl;

                for (i=0; i<points._points.cols(); ++i)
                {
                    for (j=0; j<points._points.rows(); ++j)
                    {
                        if (j > 0) {outdata << " , ";}
                        outdata << _points(j,i);
                    }
                    outdata << " , " << spin(0,i);
                    outdata << std::endl;
                }
                outdata.close();
            }
            else
            {
                std::cout << "Error: file could not be opened" << std::endl;
                std::cerr << "Error: file could not be opened" << std::endl;
                exit(1);
            }

            return;
        }

        Eigen::RowVectorXd targetgetDistanceSq(int idx)
        {
            return targetSq.block(0,std::get<0>(points.slice(0,idx)),1,std::get<1>(points.slice(0,idx)));
        }
        
        double targetgetLargestDistSq(int idx)
        {
            Eigen::RowVectorXd _Sq = targetSq.block(0,std::get<0>(points.slice(0,idx)),1,std::get<1>(points.slice(0,idx)));
            
            Eigen::RowVectorXd::Index maxIndex;
            double _ = _Sq.maxCoeff(&maxIndex);
            while (_ == std::numeric_limits<double>::infinity())
            {
                _Sq(0,maxIndex) = 0;
                _ = _Sq.maxCoeff(&maxIndex);

                if (_ == 0)
                    return std::numeric_limits<double>::infinity();
            }
            return _Sq(0,maxIndex);
        }

        double targetgetSmallestDistSq(int idx)
        {
            
            Eigen::RowVectorXd _Sq = targetSq.block(0,std::get<0>(points.slice(0,idx)),1,std::get<1>(points.slice(0,idx)));
                

            Eigen::RowVectorXd::Index minIndex;
            double _ = _Sq.minCoeff(&minIndex);
            while (_ == 0)
            {
                _Sq(0,minIndex) = std::numeric_limits<double>::infinity();
                _ = _Sq.minCoeff(&minIndex);

                if (_ == std::numeric_limits<double>::infinity())
                    return std::numeric_limits<double>::infinity();
            }
            return _Sq(0,minIndex);
        }



    };

        //----------------------------------------------------------------------------------------------------------
}