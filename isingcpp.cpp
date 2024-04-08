#include "Headers/IsingSystem.h"
#include "Headers/Tools.h"
#include <windows.h>

int main()
{
    //for (int eq : {1,2,3,4})
    for (int eq : {0})
{
    //std::vector<int> Ns = {4,8,16,32,64,128,256,512,1024,2048,4096};
    std::vector<int> Ns = {32};
    for (int N : Ns)
    {
        std::vector<double> Spins = {-1,1};
        std::vector<double> SpinsDist = {0.5,0.5};
        std::vector<double> Ts = Tools::linspace(9.5,3.5,13);
        if (Ts[0] == 0) Ts[0] = 0.01;
        if (Ts.back() == 0) Ts.back() = 0.01;
        Eigen::MatrixXd p = Tools::PointArray2D(N,N);
        Eigen::VectorXd Ls = Eigen::VectorXd::Constant(p.rows(),1);
        Eigen::MatrixXd brav = Tools::BravaisLattice::tP2D();

        bool grid = false;
        bool getStd = true;
        bool getAc = true;

        if (!grid)
        {
            Eigen::VectorXd delta1 = p.rowwise().maxCoeff() - p.rowwise().minCoeff();
            for (int i = 0; i < p.cols(); i++)
            {
                p.col(i) = p.col(i) + Eigen::Vector2d::Random()/4;
                for (int j=0; j < p.rows(); j++)
                {
                    if (p(j,i) >= N) {p(j,i) -= N;}
                    if (p(j,i) <  0) {p(j,i) += N;}
                }
            }
            Eigen::VectorXd delta2 = p.rowwise().maxCoeff() - p.rowwise().minCoeff();
            Eigen::VectorXd delta = delta1 - delta2;
            for (int j = 0; j < p.rows(); j++)
            {
                Ls(j) = std::max(0.,Ls(j) + delta(j));
            }
        }

        Ising::Ising newObj = Ising::Ising(p, 1, 1.5, true, Ls, Spins, SpinsDist, 2,brav,eq>=3,0.,getStd,getAc);

        if ((eq == 2) || (eq == 4)) {newObj.spin = (newObj.spin.array() == newObj.spin.array()).select(Spins[0],newObj.spin);}

        std::vector<std::tuple<int,int,bool,bool,bool,std::string,bool>> Runs;

        //       n_runs, n_runs_abs, store_data, store_time, store_state, "name", update system
        Runs = {
              //{ 0   , 1         , false     , false     , true       , "ImStart"       , false }
            //{ 0     , 100       , false     , true      , false      , "UnEqTime"       , false}
            { 32   , 0         , false     , false     , true       , "Recorded"       , true }
              //,{ 1024  ,  0        , false     , false     , false       , "Initial"       , true }
              //,{ 2048  , 0         , true      , false     , false      , "Equilibriated"  , true }
        //      ,{ 0     , 100       , false     , true      , false      , "EqTime"         , false}
             // ,{ 0   , 1         , false     , false     , true       , "ImEnd"       , false }
               };
        

        std::string StoreDir;
        std::string gridStr = "ON";
        if (!grid) {gridStr = "OFF";}
        StoreDir = "LRIM"+gridStr + std::to_string(N)+"/";
        if (eq > 0) {StoreDir = "LRIM"+gridStr + "Timed"+"/";}

       {
        std::string _stemp = StoreDir.substr(0,StoreDir.size()-1);
        std::wstring stemp = std::wstring(_stemp.begin(), _stemp.end());
        LPCWSTR sw = stemp.c_str();
        CreateDirectoryW(sw,NULL);
       }
       
        if (eq > 0) {StoreDir = StoreDir+"VariedOut/";}
        
        newObj.StoreStateOut = StoreDir+"States/IsingState";
        newObj.StoreStateEveryN = 1;
        newObj.StoreDataEveryN = 4;
        newObj.StoreTimeEveryN = 1;
        
       {
        std::string _stemp = StoreDir.substr(0,StoreDir.size()-1);
        std::wstring stemp = std::wstring(_stemp.begin(), _stemp.end());
        LPCWSTR sw = stemp.c_str();
        CreateDirectoryW(sw,NULL);
       }
       {
        std::string _stemp = StoreDir+"States";
        std::wstring stemp = std::wstring(_stemp.begin(), _stemp.end());
        LPCWSTR sw = stemp.c_str();
        CreateDirectoryW(sw,NULL);
       }
       {
        std::string _stemp = StoreDir+"Images";
        std::wstring stemp = std::wstring(_stemp.begin(), _stemp.end());
        LPCWSTR sw = stemp.c_str();
        CreateDirectoryW(sw,NULL);
       }

        newObj.StoreData = false;
        newObj.StoreTime = false;
        newObj.StoreState = false;
        if (eq == 0) {newObj.RunCyclesOver(0,Ts[0],true);}
        newObj.CycleCount = 0;
        for (double T : Ts)
        {
            for (std::tuple<int,int,bool,bool,bool,std::string,bool> info : Runs)
            {
                if (std::get<0>(info) + std::get<1>(info) == 0) {continue;}
                std::cout << "Temperature " << T << ": Running  " + std::get<5>(info) << '\n';
                newObj.StoreData = std::get<2>(info);
                newObj.StoreTime = std::get<3>(info);
                newObj.StoreState = std::get<4>(info);

                newObj.RunCyclesOver(std::get<0>(info),T,std::get<6>(info));
                newObj.RunCyclesFor(std::get<1>(info),T,std::get<6>(info));
                
                std::cout << "Temperature " << T << ": Finished " + std::get<5>(info) << '\n';
            }
            newObj.Reset(false);
            std::cout << "Temperature " << T << " Completed" << '\n';
        }

        if (newObj.M.size() > 0)
        {

            std::string headers = "T , M , E , X , C, |M| , X(|M|) , Run";
            if (getStd) {headers = "T , M , E , X , C, |M| , X(|M|) , MErr , EErr , XErr , CErr , |M|Err , X(|M|)Err , Run";}
            if (getStd && getAc) {headers = "T , M , E , X , C, |M| , X(|M|) , MErr , EErr , XErr , CErr , |M|Err , X(|M|)Err , MAC , EAC , XAC , CAC , |M|AC , X(|M|)AC , Run";}
            std::vector<long double> _Ts;
            std::vector<std::string> _Run;
            
            for (double T : Ts)
            {
                for (std::tuple<int,int,bool,bool,bool,std::string,bool> info : Runs)
                {
                    if (std::get<0>(info) + std::get<1>(info) == 0) {continue;}
                    if (std::get<2>(info))
                    {
                        _Ts.push_back(T);
                        _Run.push_back(std::get<5>(info));
                    }
                }
            }
            int outsize = 7;
            if (getStd) {outsize = 13;}
            if (getStd && getAc) {outsize = 19;}
            
            Tools::DataMat out(outsize,_Ts.size());
            Tools::DataMatString outstr(1,_Ts.size());

            out.row(0)= Eigen::Map<Tools::DataRowVector>(_Ts.data(),1,_Ts.size());

            out.row(1)= Eigen::Map<Tools::DataRowVector>(newObj.M.data(),1,_Ts.size());
            out.row(2)= Eigen::Map<Tools::DataRowVector>(newObj.E.data(),1,_Ts.size());
            out.row(3)= Eigen::Map<Tools::DataRowVector>(newObj.X.data(),1,_Ts.size());
            out.row(4)= Eigen::Map<Tools::DataRowVector>(newObj.C.data(),1,_Ts.size());
            out.row(5)= Eigen::Map<Tools::DataRowVector>(newObj.Ma.data(),1,_Ts.size());
            out.row(6)= Eigen::Map<Tools::DataRowVector>(newObj.Xa.data(),1,_Ts.size());
            if (getStd) {
            out.row(7)= Eigen::Map<Tools::DataRowVector>(newObj.Me.data(),1,_Ts.size());
            out.row(8)= Eigen::Map<Tools::DataRowVector>(newObj.Ee.data(),1,_Ts.size());
            out.row(9)= Eigen::Map<Tools::DataRowVector>(newObj.Xe.data(),1,_Ts.size());
            out.row(10)= Eigen::Map<Tools::DataRowVector>(newObj.Ce.data(),1,_Ts.size());
            out.row(11)= Eigen::Map<Tools::DataRowVector>(newObj.Mae.data(),1,_Ts.size());
            out.row(12)= Eigen::Map<Tools::DataRowVector>(newObj.Xae.data(),1,_Ts.size());
            }
            if (getStd && getAc) {
            out.row(13)= Eigen::Map<Tools::DataRowVector>(newObj.MAC.data(),1,_Ts.size());
            out.row(14)= Eigen::Map<Tools::DataRowVector>(newObj.EAC.data(),1,_Ts.size());
            out.row(15)= Eigen::Map<Tools::DataRowVector>(newObj.XAC.data(),1,_Ts.size());
            out.row(16)= Eigen::Map<Tools::DataRowVector>(newObj.CAC.data(),1,_Ts.size());
            out.row(17)= Eigen::Map<Tools::DataRowVector>(newObj.MaAC.data(),1,_Ts.size());
            out.row(18)= Eigen::Map<Tools::DataRowVector>(newObj.XaAC.data(),1,_Ts.size());
            }

            outstr.row(0)= Eigen::Map<Tools::DataRowVectorString>(_Run.data(),1,_Ts.size());

            std::cout << out << '\n';
            std::cout << outstr << '\n';

            std::string ModStr = "";
            if (eq == 1) {ModStr = "Modified/";}
            if (eq == 2) {ModStr = "ModifiedOrdered/";}
            if (eq == 3) {ModStr = "BruteForce/";}
            if (eq == 4) {ModStr = "BruteForceOrdered/";}

           {
            std::string _stemp = StoreDir + ModStr.substr(0,ModStr.size()-1);
            std::wstring stemp = std::wstring(_stemp.begin(), _stemp.end());
            LPCWSTR sw = stemp.c_str();
            CreateDirectoryW(sw,NULL);
           }

            Tools::WriteMatrix(ModStr,StoreDir+ModStr+std::to_string(N)+"_cppOut.csv", headers, out, outstr);
           headers = "R";

            for (int i = 0; i < (int)_Ts.size(); i++)
            {
                headers = headers + " , ";
                headers = headers + std::to_string(_Ts[i]) + _Run[i];
            }
            if (getStd) 
                for (int i = 0; i < (int)_Ts.size(); i++)
                {
                    headers = headers + " , ";
                    headers = headers + std::to_string(_Ts[i]) + _Run[i]+"Err";
                }
            if (getStd && getAc) 
                for (int i = 0; i < (int)_Ts.size(); i++)
                {
                    headers = headers + " , ";
                    headers = headers + std::to_string(_Ts[i]) + _Run[i]+"AC";
                }

            outsize = 1+_Ts.size();
            if (getStd) {outsize = 1+2*_Ts.size();}
            if (getStd && getAc) {outsize = 1+3*_Ts.size();}
            out = Tools::DataMat(newObj.gBins.size(),outsize);

            out.col(0)= Eigen::Map<Tools::DataRowVector>(newObj.gBins.data(),1,newObj.gBins.size());

            for (int i = 0; i < (int) _Ts.size(); i++)
            {
                out.col(i+1)= Eigen::Map<Tools::DataRowVector>(newObj.g[i].data(),1,newObj.gBins.size());
                if (getStd) out.col(i+1+_Ts.size()) = Eigen::Map<Tools::DataRowVector>(newObj.ge[i].data(),1,newObj.gBins.size());
                if (getStd && getAc) out.col(i+1+2*_Ts.size()) = Eigen::Map<Tools::DataRowVector>(newObj.gAC[i].data(),1,newObj.gBins.size());
            }

            out.transposeInPlace();
            std::cout << out << '\n';
            Tools::WriteMatrix(ModStr,StoreDir+ModStr+std::to_string(N)+"_cppOutG.csv", headers, out);
        }

        if (newObj.Times.size() > 0)
        {
            std::string headers = "T , Times , Delta_Times, Run";
            std::vector<long double> _Ts;
            std::vector<std::string> _Run;
            
            for (double T : Ts)
            {
                for (std::tuple<int,int,bool,bool,bool,std::string,bool> info : Runs)
                {
                    if (std::get<0>(info) + std::get<1>(info) == 0) {continue;}
                    if (std::get<3>(info))
                    {
                        _Ts.push_back(T);
                        _Run.push_back(std::get<5>(info));
                    }
                }
            }

            Tools::DataMat out(3,_Ts.size());
            Tools::DataMatString outstr(1,_Ts.size());

            out.row(0)= Eigen::Map<Tools::DataRowVector>(_Ts.data(),1,_Ts.size());
            out.row(1)= Eigen::Map<Tools::DataRowVector>(newObj.Times.data(),1,_Ts.size());
            out.row(2)= Eigen::Map<Tools::DataRowVector>(newObj.Times_std.data(),1,_Ts.size());

            outstr.row(0)= Eigen::Map<Tools::DataRowVectorString>(_Run.data(),1,_Ts.size());

            std::cout << out << '\n';
            std::cout << outstr << '\n';

            std::string ModStr = "";
            if (eq == 1) {ModStr = StoreDir +"Modified/";}
            if (eq == 2) {ModStr = StoreDir +"ModifiedOrdered/";}
            if (eq == 3) {ModStr = StoreDir +"BruteForce/";}
            if (eq == 4) {ModStr = StoreDir +"BruteForceOrdered/";}

          {
            std::string _stemp = ModStr.substr(0,ModStr.size()-1);
            std::wstring stemp = std::wstring(_stemp.begin(), _stemp.end());
            LPCWSTR sw = stemp.c_str();
            CreateDirectoryW(sw,NULL);
          }
            Tools::WriteMatrix(ModStr,std::to_string(N)+"_cppTimerOut.csv", headers, out, outstr);
        }
    }
    
}
	return 0;
}

