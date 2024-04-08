
#include "PointsSplitter.h"
#include <future>
#include <vector>
#include <iostream>
#include <numeric>

using Eigen::all;
using Eigen::last;
using PointsSplitter::SplitterObject;

namespace MCSys
{

        //----------------------------------------------------------------------------------------------------------

    class MCSystem
    {
    public:

        //----------------------------------------------------------------------------------------------------------

        SplitterObject points;
        int pointIdx;
        double iT;
        bool UpdateSystem;

        int CycleCount;
        bool StoreData;
        bool StoreState;
        bool StoreTime;
        std::string StoreStateOut;
        int StoreStateEveryN = 1;
        int StoreDataEveryN = 1;
        int StoreTimeEveryN = 1;

        //----------------------------------------------------------------------------------------------------------
        
        unsigned long long last_time_recorded;
        std::vector<long double> TimeStore;
        std::vector<long double> Times;
        std::vector<long double> Times_std;
        bool LongCycle;

        //----------------------------------------------------------------------------------------------------------

        Eigen::Vector2d dE;
        Eigen::Matrix2Xd energyBounds;
        Eigen::RowVectorXi depths;
        
        //----------------------------------------------------------------------------------------------------------

        std::vector<std::function<Eigen::Vector2d (int)>> DepthFunctions;
        int max_depth;
        bool SplitInitially;

        //----------------------------------------------------------------------------------------------------------

        MCSystem()
        {
            return;
        }

        void Reset(bool systemReset = true)
        {
            ResetBody(systemReset);

            TimeStore = {};
            resetTime();

            CycleCount = 0;
        }
        

        void BuildPoints(Eigen::MatrixXd _points,bool _PBCImaging,Eigen::Vector2d _PBCSpacing,int recDepth,int splitCount, Eigen::MatrixXd tMatrix)
        {
            srand(time(NULL));
            points = SplitterObject(_points, _PBCImaging, _PBCSpacing, recDepth, splitCount, tMatrix);
            Reset();
            /*
            PreCalculate();
            Reset();
            points.ResetMats();
            */
        }

        //----------------------------------------------------------------------------------------------------------

        void RunCyclesFor(int count,double T, bool _UpdateSystem)
        {
            UpdateSystem = _UpdateSystem;
            LongCycle = false;
            if (count < 1) {return;}
            int _j = 0;
            while (_j < count){
                if (StoreState && _j % StoreStateEveryN == 0) {WriteState(StoreStateOut,T);}
                if (StoreData && _j % StoreDataEveryN == 0) {StoreInfo();}
                if (StoreTime && _j % StoreTimeEveryN == 0) {StoreTimes();}
                CycleFor(1,T);
                _j++;
            }
            if (StoreTime) {SaveTime();}
            if (StoreData) {SaveInfo(T);}
        }

        void RunCyclesOver(int count,double T, bool _UpdateSystem)
        {
            UpdateSystem = _UpdateSystem;
            LongCycle = true;
            if (count < 1) {return;}
            int _j = 0;
            while (_j < count){
                if (StoreState && _j % StoreStateEveryN == 0) {WriteState(StoreStateOut,T);}
                if (StoreData && _j % StoreDataEveryN == 0) {StoreInfo();}
                if (StoreTime && _j % StoreTimeEveryN == 0) {StoreTimes();}
                CycleFor(points._points.cols(),T);
                _j++;
            }
            if (StoreTime) {SaveTime();}
            if (StoreData) {SaveInfo(T);}
        }

        void CycleFor(int count, double T)
        {
            if (count < 1) {return;}
            iT = 1/T;
            int __i = 0;
            while (__i < count)
            {
                double p = -(std::log(static_cast<double>(rand()) / static_cast<double>(RAND_MAX)) / iT);
                int _pointIdx = -1;
                if (LongCycle)
                    _pointIdx = __i;
                else
                    while (_pointIdx >= points._points.cols() || _pointIdx < 0 )
                        _pointIdx = (int) round((static_cast<double>(rand()) / static_cast<double>(RAND_MAX)) * (points._points.cols()+2))-1;
                CollapseFor(_pointIdx, p);
                if (UpdateSystem) {CycleCount++;}
                __i++;
            }
        }

        //----------------------------------------------------------------------------------------------------------

        void CollapseFor(int _pointIdx, double p)
        {
            pointIdx = _pointIdx;
            points.Initialise(pointIdx);
            Initialise();
            InitHandle();
            Eigen::Vector2d eb = TotalEnergyBounds();
            while (true)
            {
                if (p >= eb(1, 0))
                {
                    if (UpdateSystem) {OnSuccess();}
                    return;
                }
                else if (p < eb(0, 0))
                {
                    if (UpdateSystem) {OnFail();}
                    return;
                }
                else
                {
                    eb = NextEnergyBounds();
                }
            }
        }

        //----------------------------------------------------------------------------------------------------------

        void InitHandle()
        {
            energyBounds = Eigen::Vector2d::Constant(0);
            depths = Eigen::RowVectorXi::Zero(1);
            if (SplitInitially)
            {
                SplitSection(0);
                depths  = Eigen::RowVectorXi::Zero(points.dimConst);
                InitialiseSection(0);
            }
            else
                CalcEnergyBounds(0);
        }
        
        //----------------------------------------------------------------------------------------------------------

        Eigen::Vector2d TotalEnergyBounds()
        {
            dE = energyBounds.rowwise().sum();
            return dE;
        }

        //----------------------------------------------------------------------------------------------------------

        Eigen::Vector2d NextEnergyBounds()
        {
            SplitNext();
            return TotalEnergyBounds();
        }

        void SplitNext()
        {
            Eigen::MatrixXd::Index maxIndex;
            (energyBounds.row(1) - energyBounds.row(0)).maxCoeff(&maxIndex);
            if (depths(0,maxIndex) == (int)DepthFunctions.size())
            {
                depths(0,maxIndex) = depths(0,maxIndex)+1;
                CalcEnergyBounds(maxIndex);
            }
            else
            {
                SplitSection(maxIndex);
                InitialiseSection(maxIndex);
            }
        }

        void UpdateHandleMatricesSplit(int idx)
        {
            energyBounds = ExtendMatDouble(idx, energyBounds, 0,points.dimConst);
            depths = ExtendMatInt(idx, depths, depths(0, idx) + 1,points.dimConst);
        }

        void SplitSection(int idx)
        {
            UpdateSplitMatrices(idx);
            UpdateHandleMatricesSplit(idx);
            points.splitOnIndex(idx);
        }
        void InitialiseSection(int idx)
        {
            /*
            for (int x = idx; x < idx+ (int) points.dimConst; x++)
                CalcEnergyBounds(x);
            */
            std::vector<std::future<void>> futs;
            for (int x = idx; x < idx+ (int) points.dimConst; x++)
                futs.push_back(std::async(std::launch::async, CalcEnergyBounds, this, x));
            for (int x = 0; x < (int) points.dimConst; x++)
                futs[x].get();
        }

        //----------------------------------------------------------------------------------------------------------

        void CalcEnergyBounds(int idx)
        {
            int s = std::get<1>(points.slice(0,idx));
            if (s == 0)
            {
                energyBounds(all, idx) = Eigen::Vector2d::Constant(0);
            }
            else if (s == 1 && std::get<0>(points.slice(0,idx)) == pointIdx)
            {
                energyBounds(all, idx) = Eigen::Vector2d::Constant(0);
            }
            else 
            {
                energyBounds(all, idx) = DepthFunctions[depths(0, idx)](idx);
            }
            return;
        }

        //----------------------------------------------------------------------------------------------------------

        Eigen::MatrixXd ExtendMatDouble(int idx, Eigen::MatrixXd Mat, double _toIns, int n)
        {
            Eigen::MatrixXd _Mat(Mat.rows(), Mat.cols() + n -1);
            if (idx > 0)
            {
                _Mat.leftCols(idx) = Mat.leftCols(idx);
            }
            _Mat(all, Eigen::seq(idx, idx + n - 1)) = Eigen::MatrixXd::Constant(Mat.rows(), n, _toIns);
            if (Mat.cols() - 1 > idx)
            {
                _Mat.rightCols(Mat.cols()-idx-1) = Mat.rightCols(Mat.cols()-idx-1);
            }
            return _Mat;
        }

        Eigen::MatrixXi ExtendMatInt(int idx, Eigen::MatrixXi Mat, int _toIns, int n)
        {
            Eigen::MatrixXi _Mat(Mat.rows(), Mat.cols() + n - 1);
            if (idx > 0)
            {
                _Mat.leftCols(idx) = Mat.leftCols(idx);
            }
            _Mat(all, Eigen::seq(idx, idx + n - 1)) = Eigen::MatrixXi::Constant(Mat.rows(), n, _toIns);
            if (Mat.cols() - 1 > idx)
            {
                _Mat.rightCols(Mat.cols()-idx-1) = Mat.rightCols(Mat.cols()-idx-1);
            }
            return _Mat;
        }

        //----------------------------------------------------------------------------------------------------------

        void SaveTime()
        {
            long double T = mean(TimeStore);
            long double Tstd = StDev(TimeStore,T);
            Times.push_back(T);
            Times_std.push_back(Tstd);

            TimeStore = {};
            resetTime();
        }

        void StoreTimes()
        {
            TimeStore.push_back(getTimePassed());
        }
        
        long double getTimePassed()
        {
            unsigned long long t;
            unsigned long long _t = last_time_recorded;
            t = std::chrono::system_clock::now().time_since_epoch() / std::chrono::microseconds(1);
            last_time_recorded = t;
            return (long double) (t - _t)/1000000;
        }

        void resetTime()
        {
            last_time_recorded = std::chrono::system_clock::now().time_since_epoch() / std::chrono::microseconds(1);
        }

        //----------------------------------------------------------------------------------------------------------

        virtual void UpdateSplitMatrices(int idx) { return; }
        virtual void Initialise() { return; }
        virtual void OnSuccess() { return; }
        virtual void OnFail() { return; }
        virtual void StoreInfo() { return; }
        virtual void SaveInfo(double T) { return; }
        virtual void ResetBody(bool systemReset = true) { return; }
        virtual void PreCalculate() { return; }
        virtual void WriteState(std::string file, double T) { return; }

        //----------------------------------------------------------------------------------------------------------

        std::vector<long double> prod(std::vector<long double> vector1,std::vector<long double> vector2)
        {
            std::vector<long double> accum = {};
            for (int i = 0; i < (int)vector1.size(); i++)
                accum.push_back(vector1[i]*vector2[i]);
            return accum;
        }

        long double mean(std::vector<long double> vector)
        {
            long double sum = std::accumulate(std::begin(vector), std::end(vector), 0.0);
            long double mean =  sum / vector.size();
            return mean;
        }

        long double StDev(std::vector<long double> vector, long double mean)
        {
            long double accum = 0.0;
            std::for_each (std::begin(vector), std::end(vector), [&](const double d) {
                accum += (d - mean) * (d - mean);
            });
            long double stdev = sqrt(accum / ((long double)(vector.size()-1)));
            return stdev;
        }
        
        std::vector<long double> JackKnife(std::function<long double (std::vector<long double>,std::vector<long double>)> fun,std::vector<long double> V1,std::vector<long double> V2)
        {
            std::vector<long double> accum = {};
            for (int i = 0; i < (int) V1.size(); i++)
            {
                std::vector<long double> V1_ = {};
                std::vector<long double> V2_ = {};
                
                for (int j = 0; j < (int) V1.size(); j++)
                {
                    if (!(j==i))
                    {
                        V1_.push_back(V1[j]);
                        V2_.push_back(V2[j]);
                    }
                }
                accum.push_back(fun(V1_,V2_));
            }
            return accum;
        }

        std::vector<long double> JackKnife3(std::function<long double (std::vector<long double>,std::vector<long double>,std::vector<long double>)> fun,std::vector<long double> V1,std::vector<long double> V2,std::vector<long double> V3)
        {
            std::vector<long double> accum = {};
            for (int i = 0; i < (int) V1.size(); i++)
            {
                std::vector<long double> V1_ = {};
                std::vector<long double> V2_ = {};
                std::vector<long double> V3_ = {};
                
                for (int j = 0; j < (int) V1.size(); j++)
                {
                    if (!(j==i))
                    {
                        V1_.push_back(V1[j]);
                        V2_.push_back(V2[j]);
                        V3_.push_back(V3[j]);
                    }
                }
                accum.push_back(fun(V1_,V2_,V3_));
            }
            return accum;
        }


        long double autocorrelation(std::vector<long double> vector, long double mean,long double stddev,int lag = 1)
        {
            long double accum = 0.0;
            for (int i = 0; i+lag < (int)vector.size(); i++)
                accum += (vector[i]-mean)*(vector[i+lag]-mean);
            long double AC = accum * (1 / ((long double) (vector.size()-lag) * stddev*stddev));
            return AC;
        }
    };

        //----------------------------------------------------------------------------------------------------------

}