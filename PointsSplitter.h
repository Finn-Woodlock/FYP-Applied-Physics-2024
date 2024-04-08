#include <cmath>
#include "../ThirdParty/eigen/Eigen/Eigen"
#include <tuple>
#include <iostream>

using Eigen::all;
using Eigen::last;

typedef Eigen::Matrix<std::tuple<int,int>, 1, Eigen::Dynamic> SliceMatrix;

namespace PointsSplitter
{

        //----------------------------------------------------------------------------------------------------------
    
    class SplitterObject
    {
    public:

        //----------------------------------------------------------------------------------------------------------

        Eigen::MatrixXd _points;
        Eigen::VectorXd point;
        unsigned int dimConst;

        //----------------------------------------------------------------------------------------------------------

        SliceMatrix slice;
        Eigen::RowVectorXi SliceAccumulationIndex;
        Eigen::RowVectorXi SliceIndex;
        int RecursionColCount;
        int RecursionCount;

        //----------------------------------------------------------------------------------------------------------

        Eigen::MatrixXd T;

        //----------------------------------------------------------------------------------------------------------

        Eigen::VectorXd PBCBarrier;
        Eigen::MatrixXd T_PBCBarrier;
        Eigen::VectorXd T_PBCBarrier2;
        bool PBCImaging;

        Eigen::MatrixXd tMatrix;
        bool checkMatrix;
        Eigen::RowVectorXd Sq;

        //----------------------------------------------------------------------------------------------------------



        //----------------------------------------------------------------------------------------------------------

        SplitterObject()
        {
            return;
        }

        SplitterObject(Eigen::MatrixXd __points, bool _PBCImaging, Eigen::VectorXd _PBCSpacing, int recursionDepth, int splitCount, Eigen::MatrixXd _tMatrix)
        {
            PBCImaging = _PBCImaging;
            if (!PBCImaging){_PBCSpacing = Eigen::VectorXd::Constant(__points.rows(),0.01);}

            tMatrix=_tMatrix;
            checkMatrix = false;
            for (int i=0; i < tMatrix.cols(); i++)
            {
                for (int j=0; j < tMatrix.rows(); j++)
                {
                    if ((i==j && !(tMatrix(i,j) == 1.)) || (!(i==j) && !(tMatrix(i,j) == 0.))) {checkMatrix = true;}
                }
            }
            
            PBCBarrier = (_PBCSpacing+__points.rowwise().maxCoeff());

            _points = tMatrix*SortPoints(__points, recursionDepth, splitCount);
            
            if (checkMatrix)
            {
                T_PBCBarrier = Eigen::MatrixXd((int) (std::pow(3,_points.rows())-1)/2,_points.rows());
                T_PBCBarrier2 = Eigen::VectorXd((int) (std::pow(3,_points.rows())-1)/2);
                int row = 0;
                for (int i = 0; i < _points.rows(); i++)
                {
                    for (int j = 0; j < std::pow(3,_points.rows()-1-i); j++)
                    {
                        Eigen::VectorXd mPBCBarrier = Eigen::VectorXd::Zero(_points.rows());
                        mPBCBarrier(i,0) = PBCBarrier(i,0);
                        for (int k = 1; k+i < _points.rows();k++)
                        {
                            mPBCBarrier(k+i,0) = (((int) std::floor(((double)j)/((double)std::pow(3,k-1)))%3)-1) * PBCBarrier(k+i,0);
                        }
                        T_PBCBarrier.row(row) = (tMatrix * mPBCBarrier).transpose();
                        T_PBCBarrier2(row,0) = T_PBCBarrier.row(row) * T_PBCBarrier.row(row).transpose();
                        row++;
                    }
                }
            }

            ResetMats();
        }

        //----------------------------------------------------------------------------------------------------------
        
        void Initialise(int pointIdx)
        {
            point = _points(all, pointIdx);
            //Sq = PBCImage(_points);
            //points = PBCImage(_points);
            Sq = SqPBCImage(_points);
            Sq(0,pointIdx) = std::numeric_limits<double>::infinity();
            ResetMats();
        }

        void ResetMats()
        {
            Eigen::RowVectorXi _SliceIndex(2);
            _SliceIndex << 0, RecursionCount;
            SliceIndex = _SliceIndex;
            SliceMatrix _slice(1);
            std::tuple<int,int> _slice1 = {0,SliceAccumulationIndex(0,last)};
            _slice << _slice1;
            slice = _slice;
        }

        //----------------------------------------------------------------------------------------------------------

        Eigen::MatrixXd SortPoints(Eigen::MatrixXd __points, int recursionDepth, int splitCount)
        {
            dimConst = std::pow(splitCount,__points.rows());

            T = Eigen::MatrixXd::Zero(__points.rows()*2, __points.rows()*2*dimConst);
            for (float idx = 0; idx < dimConst;idx++)
            {
                Eigen::MatrixXd _T = Eigen::MatrixXd::Zero(2*__points.rows(), 2*__points.rows());
                for (unsigned int i = 0; i < __points.rows(); i++)
                {
                    for (float j = 0; j < splitCount; j++)
                    {
                        if (((int)idx%(int)std::pow(splitCount,i+1))/(int)std::pow(splitCount,i)==(int)j)
                        {
                            _T(i, i) = (splitCount-j)/splitCount;
                            _T(i, i+__points.rows()) = j/splitCount;
                            _T(i+__points.rows(), i) = (splitCount-j-1)/splitCount;
                            _T(i+__points.rows(), i+__points.rows()) = (j+1)/splitCount;
                        }
                    }
                }
                T.block(0,idx*__points.rows()*2,2*__points.rows(),2*__points.rows()) = _T;
            }


            Eigen::VectorXd bounds1(__points.rows()*2);
            bounds1 << Eigen::VectorXd::Zero(__points.rows()), PBCBarrier;

            SliceAccumulationIndex = Eigen::RowVectorXi::Zero(std::pow(splitCount,__points.rows()*recursionDepth) + 1);
            RecursionCount=0;
            RecursionColCount=0; 
            Eigen::MatrixXd out = SortBox(__points,bounds1,recursionDepth);
            return out;
        }

        Eigen::MatrixXd SortBox(Eigen::MatrixXd __points, Eigen::VectorXd _bounds, int recursion_depth)
        {
            if (recursion_depth == 0)
            {
                RecursionColCount = RecursionColCount + __points.cols();
                SliceAccumulationIndex(0,RecursionCount+1) = RecursionColCount;
                RecursionCount++;
                return __points;
            }
            else
            {
                Eigen::MatrixXd _T(2*__points.rows(),2*__points.rows());
                Eigen::MatrixXd out(__points.rows(),__points.cols());
                int cols = 0;
                for (unsigned int idx = 0; idx < dimConst;idx++)
                {
                    _T = T.block(0,idx*__points.rows()*2,__points.rows()*2,__points.rows()*2);
                    Eigen::MatrixXd out1 = SortBox(__points(all,WhereInBounds(__points,_T * _bounds)),_T * _bounds,recursion_depth-1);
                    out(all,Eigen::seq(cols,cols+out1.cols()-1)) = out1;
                    cols+=out1.cols();
                }
                return out;
            }
        }

        std::vector<int> WhereInBounds(Eigen::MatrixXd __points,Eigen::VectorXd bound)
        {
            std::vector<int> out = {};
            for (int i = 0; i < __points.cols(); i++)
            {
                bool Success = true;
                for (int j = 0; j < __points.rows(); j++)
                    Success = Success && __points(j,i) < bound(__points.rows() + j,0) && __points(j,i) >= bound(j,0);
                if (Success)
                {
                    out.push_back(i);
                }
            }
            return out;
        }

        //----------------------------------------------------------------------------------------------------------

        Eigen::RowVectorXd SqPBCImageOverride(Eigen::MatrixXd __points, int idx)
        {
            Eigen::VectorXd m_point = point;
            point = __points(all,idx);
            Eigen::RowVectorXd _Sq = SqPBCImage(__points,true);
            point = m_point;
            return _Sq;
        }

        Eigen::RowVectorXd SqPBCImage(Eigen::MatrixXd __points, bool Reduce = true)
        {
            Eigen::RowVectorXd out(__points.cols());
            Eigen::MatrixXd Relpoints = (__points.colwise() - point);
            if (PBCImaging)
            {
                if (!checkMatrix)
                {
                    Eigen::RowVectorXd Zers = Eigen::RowVectorXd::Zero(__points.cols());
                    if (Reduce)
                    {
                        for (int row = 0; row < __points.rows(); row++)
                        {
                            Relpoints.row(row) += (Relpoints.row(row).array() < -PBCBarrier(row,0)/2).select(PBCBarrier(row,0),Zers) - (Relpoints.row(row).array() > PBCBarrier(row,0)/2).select(PBCBarrier(row,0),Zers);
                        }
                    }
                    out = Relpoints.array().square().colwise().sum();
                    
                }
                else
                {
                    out = Relpoints.array().square().colwise().sum();
                    if (Reduce)
                    {
                        Eigen::MatrixXd BestSq = Eigen::MatrixXd::Zero(2,__points.cols());
                        for (int row = 1; row-1 < T_PBCBarrier.rows(); row++)
                        {
                            Eigen::RowVectorXd PMat = T_PBCBarrier.row(row-1);
                            double PCond = T_PBCBarrier2(row-1,0);
                            BestSq.row(1) = Eigen::RowVectorXd::Constant(__points.cols(),PCond)-(2*PMat*Relpoints).array().abs().matrix();
                            BestSq.row(0) = BestSq.colwise().minCoeff();
                        }
                        out += BestSq.row(0);
                    }
                }
            }
            else
            {
                out = Relpoints.array().square().colwise().sum();
            }
            return out;
        }

        Eigen::MatrixXd PBCImage(Eigen::MatrixXd __points, bool Reduce = true)
        {
            Eigen::MatrixXd Relpoints = (__points.colwise() - point);
            if (PBCImaging)
            {
                if (!checkMatrix)
                {
                    Eigen::RowVectorXd Zers = Eigen::RowVectorXd::Zero(__points.cols());
                    if (Reduce)
                    {
                        for (int row = 0; row < __points.rows(); row++)
                        {
                            Relpoints.row(row) += (Relpoints.row(row).array() < -PBCBarrier(row,0)/2).select(PBCBarrier(row,0),Zers) - (Relpoints.row(row).array() > PBCBarrier(row,0)/2).select(PBCBarrier(row,0),Zers);
                        }
                    }
                    
                }
                else
                {
                    Eigen::RowVectorXd Zers = Eigen::RowVectorXd::Zero(__points.cols());
                    if (Reduce)
                    {
                        Eigen::MatrixXd BestSq = Eigen::MatrixXd::Zero(3,__points.cols());
                        Eigen::RowVectorXi BestX = Eigen::RowVectorXi::Zero(__points.cols());
                        for (int row = 1; row-1 < T_PBCBarrier.rows(); row++)
                        {
                            Eigen::RowVectorXd PMat = T_PBCBarrier.row(row-1);
                            double PCond = T_PBCBarrier2(row-1,0);
                            BestSq.row(1) = Eigen::RowVectorXd::Constant(__points.cols(),PCond)+2*PMat*Relpoints;
                            BestSq.row(2) = Eigen::RowVectorXd::Constant(__points.cols(),PCond)-2*PMat*Relpoints;
                            BestX = (BestSq.row(0).array() < BestSq.row(1).array()).select(BestX,row);
                            BestX = (BestSq.row(0).array() < BestSq.row(2).array()).select(BestX,-row);
                            BestSq.row(0) = BestSq.colwise().minCoeff();
                        }
                        for (int row = 1; row-1 < T_PBCBarrier.rows(); row++)
                        {
                            Eigen::RowVectorXd PMat = T_PBCBarrier.row(row-1);
                            for (int i = 0; i < __points.rows(); i++)
                            {
                                __points.row(i) += (BestX.array()==row).select(PMat(0,i),Zers);
                                __points.row(i) += (BestX.array()==-row).select(-PMat(0,i),Zers);
                            }
                        }
                    }
                }
            }
            return __points;
        }
        
        //----------------------------------------------------------------------------------------------------------

        Eigen::RowVectorXd getDistanceSq(int idx)
        {
            return Sq.block(0,std::get<0>(slice(0,idx)),1,std::get<1>(slice(0,idx)));
        }
        
        double getLargestDistSq(int idx)
        {
            Eigen::RowVectorXd _Sq = Sq.block(0,std::get<0>(slice(0,idx)),1,std::get<1>(slice(0,idx)));
            
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

        double getSmallestDistSq(int idx)
        {
            
            Eigen::RowVectorXd _Sq = Sq.block(0,std::get<0>(slice(0,idx)),1,std::get<1>(slice(0,idx)));
                

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

        //----------------------------------------------------------------------------------------------------------

        void splitOnIndex(int idx)
        {
            Eigen::RowVectorXi _SliceIndex(SliceIndex.cols()+dimConst-1);
            SliceMatrix _slice(slice.cols()+dimConst-1);

            int dIndex = (SliceIndex(0,idx+1) - SliceIndex(0,idx))/dimConst;

            if (idx > 0)
            {
                _SliceIndex.leftCols(idx) = SliceIndex.leftCols(idx);
                _slice.leftCols(idx) = slice.leftCols(idx);
            }

            for (unsigned int i = 0; i < dimConst; i++)
            {
                int Id = SliceIndex(0,idx)+i*dIndex;
                _SliceIndex(0, idx+i) = Id;
                std::tuple<int,int> sl = {SliceAccumulationIndex(0,Id),
                                          SliceAccumulationIndex(0,Id+dIndex)-SliceAccumulationIndex(0,Id)};
                _slice(0, idx+i) = sl;
            }

            if (slice.cols() - 1 > idx)
            {
                _slice.rightCols(slice.cols()-idx-1) = slice.rightCols(slice.cols()-idx-1);
            }
            _SliceIndex.rightCols(SliceIndex.cols()-idx-1) = SliceIndex.rightCols(SliceIndex.cols()-idx-1);

            SliceIndex = _SliceIndex;
            slice = _slice;
            return;
        }

        //----------------------------------------------------------------------------------------------------------

    };

        //----------------------------------------------------------------------------------------------------------

}