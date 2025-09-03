#undef VSSFOAM_DEBUG

#include <iomanip>

#include "soilGlobal.h"
#include "soilMath.h"
#include "fvCFD.H"
#include "pisoControl.H"
#include "argList.H"

#include "wallFvPatch.H"
#include "OFstream.H"

#ifdef ADAPTIVE_MESH
#include "dynamicFvMesh.H"
#endif

#include "retentionDataInterface.h"
#include "retentionDataVg.h"
#include "retentionDataKs.h"
#include "retentionDataBc.h"
#include "retentionDataExp.h"
#include "retentionDataHvk.h"
#include "retentionDataFilmVg.h"
#include "retentionDataFilmKs.h"
#include "retentionDataFilmVaporNsVg.h"
#include "retentionDataFilmVaporNsKs.h"
#include "retentionDataFilmVaporStdVg.h"
#include "retentionDataVaporStd.h"
#include "simConfig.h"
#include "solverFunctions.h"

using namespace Soil::RetentionModels;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define get_str(str) #str

int main(int argc, char *argv[])
{  
    argList::addBoolOption("shp-debug", "Save SHP debug information and exit", false);
    argList::addBoolOption("skip-mb-check", "Skip mass-balance check", false);
    argList::addBoolOption("force-flushing-output", "If enabled the OFstream flushing is forced at each timestep", false);

#include "setRootCase.H"
#include "createTime.H"

#ifdef ADAPTIVE_MESH
#include "createDynamicFvMesh.H"
#else
#include "createMesh.H"
#endif

    IOdictionary transportProperties(
        IOobject(
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE)
    );

    simConfig config(transportProperties);    
    const bool is_shp_debug = args.found("shp-debug");
    const bool is_skip_mbs = args.found("skip-mb-check");
    const bool is_force_flush = args.found("force-flushing-output");
    config.algorithm.setShpDebug(is_shp_debug);
    




    List<List<scalar>> deltaTSchedule;
    runTime.controlDict().readIfPresent("deltaTSchedule", deltaTSchedule);
    int deltaTScheduleSize = deltaTSchedule.size();
    int activeDeltaTSchedule;
    int activeCurrTimeThreshold;
    if (deltaTScheduleSize > 0)
    {
        activeDeltaTSchedule = 0;
        runTime.setDeltaT(deltaTSchedule[activeDeltaTSchedule][1]);
        activeCurrTimeThreshold = int(deltaTSchedule[activeDeltaTSchedule][0]);
    };

#include "createControl.H"
#include "createFields.H"

    if (config.algorithm.isCorrectBcK()) {
        Soil::Math::setZeroGradientBoundaryCondition(K);
        K.write();
        K = retentionModel->Kh(h, K);
    }
    if (config.algorithm.isCorrectBcCv()) {
        Soil::Math::setZeroGradientBoundaryCondition(Cv);
        Cv.write();
        Cv = retentionModel->Cv(h, Cv);
    }
    if (config.algorithm.isCorrectBcTheta()) {
        Soil::Math::setZeroGradientBoundaryCondition(Theta);
        Theta.write();
        Theta = retentionModel->Theta(h, Theta);
    }




#ifdef ADAPTIVE_MESH
    magGradH.write();
#endif

    // START - MBC initialization
    scalar mbcTotalInflow = 0;
    scalar mbcInitialTWC = 0;
    scalar mbcOldTimestepTWC = 0;
    OFstream LogMBC("mass_balance_check.csv");
    
    if (!is_skip_mbs) {
        forAll(mesh.C(), celli)
        {
            mbcInitialTWC += Theta[celli] * mesh.V()[celli];
        }

        reduce(mbcInitialTWC, sumOp<scalar>());

        mbcOldTimestepTWC = mbcInitialTWC;

        if (Pstream::master())
        {
            LogMBC << "Time;mbcTimestepError;mbcTimestepErrorRelative;mbcError;mbcErrorRelative;mbcErrorCelia;mbcErrorCeliaTimestep;mbcDiffTWC;mbcDiffRelativeTWC;mbcTotalInflow;mbcDiffTWCTimestep;mbcTimestepInflow;" <<config.algorithm.getMbcHeaderInfo().c_str()<<nl;
        }
    }
     // STOP - MBC initialization

    Info << "\nCalculating soil potential distribution\n"
         << nl;

    volVectorField gradK = fvc::grad(K);

    volScalarField h_old = h;
    volScalarField K_old = K;
    volScalarField Cv_old = Cv;    
    volScalarField K_tmp = K;
    
    surfaceScalarField K_tmp_surface (
        IOobject
        (
            "K_tmp_surface",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("K_tmp_surface", dimensionSet(0, 1, -1, 0, 0, 0, 0), 0.0)
    );
    if (config.algorithm.getKInterstepEstimationDoSurfaceInterpolation() == true)
    {
        K_tmp_surface = fvc::interpolate(K_tmp);
    }
    //main loop start
    while (runTime.loop())
    {
        Info<<"Main loop"<<nl;
        if (deltaTScheduleSize > 0)  //FIXME impelent ignoring the shedule based on adjustTimeStep keyword
        {
            int currTime = runTime.value();
            if (currTime >= activeCurrTimeThreshold)
            {
                activeDeltaTSchedule++;
                if (activeDeltaTSchedule >= deltaTScheduleSize)
                {
                    scalar myDeltaT;
                    runTime.controlDict().readIfPresent("deltaT", myDeltaT);
                    runTime.setDeltaT(myDeltaT);
                    runTime.controlDict().readIfPresent("endTime", activeCurrTimeThreshold);
                    Info<<"Time SET = myDeltaT"<<nl;
                }
                else
                {
                    activeCurrTimeThreshold = int(deltaTSchedule[activeDeltaTSchedule][0]);
                    runTime.setDeltaT(deltaTSchedule[activeDeltaTSchedule][1]);
                    Info<<"Time SET = deltaTSchedule[activeDeltaTSchedule][1]"<<nl;
                }
            }
        }

        Info << "Time = " << runTime.timeName() << "(" << runTime.deltaT().value() << ")" << nl;

#ifdef ADAPTIVE_MESH
        mesh.update();
#endif
        int picard_iter = 0;
        while (picard_iter <= config.algorithm.getPicardMaxIter()) 
        {
            Info<<"\tPISO outer loop no.:"<<picard_iter<<nl;
           
                
            Info<<"\th estimated (old time), min: "<<min(h.oldTime())<<", avg: "<<average(h.oldTime())<<", max: "<<max(h.oldTime())<<endl;
            Info<<"\th estimated, min: "<<min(h)<<", avg: "<<average(h)<<", max: "<<max(h)<<endl;
            h_old = h;
            K_old = K;
            Cv_old = Cv;

            int piso_iter = 0;
            while (piso.correct())
            {
                Info<<"\t\tPISO inner loop no.:"<<piso_iter<<nl;
                while (piso.correctNonOrthogonal())
                {
                    Info<<"\t\t\tcorrectNonOrthogonal - BEFORE"<<nl;
                    
                    if (config.algorithm.getKInterstepEstimationDoSurfaceInterpolation() == true)
                    {
                        fvScalarMatrix hEqn(
                            Cv * fvm::ddt(h) - fvm::laplacian(K, h) - gradK.component(gNonZeroComp) == sources
                        );
                        hEqn.solve();
                    } else {
                        fvScalarMatrix hEqn(
                            Cv * fvm::ddt(h) - fvm::laplacian(K, h) - gradK.component(gNonZeroComp) == sources
                        );
                        hEqn.solve();
                    }
                    //    p_rghEqn.solve(mesh.solver(p_rgh.select(PISO.finalInnerIter())));

                    if (piso.finalNonOrthogonalIter())
                    {
                        Info<<"\t\t\tcorrectNonOrthogonal - FINAL"<<nl;
                    }

                    Info<<"\t\t\th estimated, min: "<<min(h)<<", avg: "<<average(h)<<", max: "<<max(h)<<endl;
                    Info<<"\t\t\tcorrectNonOrthogonal - AFTER"<<nl;
                    
                }   
            
                piso_iter++;
            }

            Info<<"\th estimated, min: "<<min(h)<<", avg: "<<average(h)<<", max: "<<max(h)<<endl;

            Theta = retentionModel->Theta(h, Theta);
            Cv = retentionModel->Cv(h, Cv);
            K = retentionModel->Kh(h, K);
                        
            if (config.algorithm.isCorrectBcK()) {
                K.correctBoundaryConditions();
            };
            if (config.algorithm.isCorrectBcCv()) {
                Cv.correctBoundaryConditions();
            };
            if (config.algorithm.isCorrectBcTheta()) {
                Theta.correctBoundaryConditions();
            };
            h.correctBoundaryConditions();
            
            switch (config.algorithm.getKInterstepEstimationMethod())
            {
                case (config.algorithm.kIntEstMethod::now):
                    K_tmp = K;
                    break;
                case (config.algorithm.kIntEstMethod::nowOldArithmeticMean):
                    K_tmp = 0.5 * (K + K.oldTime());
                    break;
                case (config.algorithm.kIntEstMethod::nowOldGeometricMean):
                    K_tmp = sqrt(K * K.oldTime());
                    break;
                case (config.algorithm.kIntEstMethod::nowOldHarmonicMean):
                    K_tmp = 2.0 * (K * K.oldTime()) / (K + K.oldTime());
                    break;
            }


            if (config.algorithm.getKInterstepEstimationDoSurfaceInterpolation() == true)
            {
                K_tmp_surface = fvc::interpolate(K_tmp);
            }

            gradK = fvc::grad(K);

            if (config.processes.isDneModelSet() == true)
            {
                h += (h - h.oldTime()) * (1.0 - Foam::exp(-runTime.deltaT() / retentionModel->getDneTau()));
            }

            if (picard_iter != 0)
            {
                scalar h_avg_residual = sqrt(average(pow((h - h_old), 2))).value();
                scalar K_avg_residual = sqrt(average(pow((K - K_old), 2))).value();
                scalar Cv_avg_residual = sqrt(average(pow((Cv - Cv_old), 2))).value();
                Info<<"\tLOOP RESIDUALS => h residual: "<<h_avg_residual<<", K residual: "<<K_avg_residual<<", Cv residual: "<<Cv_avg_residual<<nl;
                if ( (h_avg_residual <= config.algorithm.getHPicardResidual() && K_avg_residual <= config.algorithm.getKPicardResidual() && Cv_avg_residual <= config.algorithm.getCvPicardResidual()) ||
                        (h_avg_residual <= config.algorithm.getHPicardResidual() && K_avg_residual <= config.algorithm.getKPicardResidual() && config.algorithm.getCvPicardResidual() == 0.0) ||
                        (h_avg_residual <= config.algorithm.getHPicardResidual() && config.algorithm.getKPicardResidual() == 0.0                    && Cv_avg_residual <= config.algorithm.getCvPicardResidual() ) ||
                        (h_avg_residual <= config.algorithm.getHPicardResidual() && config.algorithm.getKPicardResidual() == 0.0                    && config.algorithm.getCvPicardResidual() == 0.0) )
                {
                    Info<<"\tLOOP RESIDUALS ENDING CRITERIA MET ("<<++picard_iter<<")."<<nl;
                    break;
                }
            } else {
                if (config.algorithm.getHPicardResidual() == 0.0 && config.algorithm.getKPicardResidual() == 0.0 && config.algorithm.getCvPicardResidual() == 0.0) 
                {
                    Info<<"\tLOOP RESIDUALS (Picard iteration not used)."<<nl;
                    break;
                }
            }

            if (picard_iter == config.algorithm.getPicardMaxIter())
            {
                Info<<"\tLOOP RESIDUALS NOT CONVERGED! ("<<picard_iter<<")."<<nl;
            }
                
            picard_iter++;
        }
        
        // START - MBC timestep evaluation
        
        if (!is_skip_mbs) {
            scalar mbcTimestepTWC = 0;
            scalar mbcTimestepInflow = 0;
            forAll(mesh.C(), celli)
            {
                mbcTimestepTWC += Theta[celli] * mesh.V()[celli];
            }
            reduce(mbcTimestepTWC, sumOp<scalar>());

            // evaluate boundary fluxes
            const fvPatchList &patches = mesh.boundary();
            forAll(patches, patchi)
            {
                const fvPatch &currPatch = patches[patchi];
                // Info<<"(MBC)     "<<"Patch name: "<<currPatch.name()<<nl;
                if (isType<wallFvPatch>(currPatch) == true)
                {
                    // const UList<label> &bfaceCells = mesh.boundaryMesh()[patchi].faceCells();
                    // Info<<"(MBC)     "<<"    bfaceCells.size():"<<bfaceCells.size()<<nl;
                    // Info<<"(MBC)     "<<"    bfaceCells:"<<bfaceCells<<nl;
                    const scalarField gNorm(currPatch.nf() & (g.value() / mag(g.value())));
                    const scalarField &hWall = static_cast<scalarField>(h.boundaryField()[patchi]);
                    const scalarField &hWallOld = static_cast<scalarField>(h_old.boundaryField()[patchi]);
                    const scalarField &KWall = static_cast<scalarField>(K.boundaryField()[patchi]);
                    const scalarField &KWallOld = static_cast<scalarField>(K_old.boundaryField()[patchi]);
                    forAll(currPatch, facei)
                    {
                        label faceCelli = currPatch.faceCells()[facei];
                        scalar grad = (h[faceCelli] - hWall[facei]) * currPatch.deltaCoeffs()[facei] + gNorm[facei];
                        scalar gradOld = (h_old[faceCelli] - hWallOld[facei]) * currPatch.deltaCoeffs()[facei] + gNorm[facei];
                        //scalar flux = (grad * (KWall[facei]+K[faceCelli])/2.0 + gradOld * (KWallOld[facei]+K_old[faceCelli])/2.0)/2.0;
                        //scalar flux = (grad * KWall[facei] + gradOld * KWallOld[facei])/2.0;
                        scalar flux = grad * KWall[facei];
                        scalar waterInflow = flux * runTime.deltaT().value() * currPatch.magSf()[facei];
                        #ifdef VSSFOAM_DEBUG 
                            //Info << "(MBCTIMESTEP)     " << runTime.value() << " patchi:" << patchi << "  facei:" << facei << "  grad:" << grad << " flux:" << flux << " waterInflow:" << waterInflow << " KWall:" << KWall[facei] << " K:" << K[facei] << nl;
                        #endif
                    
                        mbcTimestepInflow += waterInflow;
                    }
                }
            }
            reduce(mbcTimestepInflow, sumOp<scalar>());
            mbcTotalInflow += mbcTimestepInflow;
            reduce(mbcTotalInflow, sumOp<scalar>());
            scalar mbcError = mbcTotalInflow - ( mbcInitialTWC - mbcTimestepTWC );
            scalar mbcDiffTWC = mbcInitialTWC - mbcTimestepTWC;
            scalar mbcDiffTWCTimestep = mbcOldTimestepTWC - mbcTimestepTWC;

            scalar mbcErrorRelative = -999;
            scalar mbcDiffRelativeTWC = -999;
            if (mbcInitialTWC != 0.0) {   
                mbcErrorRelative = mbcError / mbcInitialTWC;
                mbcDiffRelativeTWC = mbcDiffTWC / mbcInitialTWC;
            }

            scalar mbcErrorCelia = -999;
            if (mbcTimestepInflow != 0.0) {
                mbcErrorCelia = (mbcInitialTWC - mbcTimestepTWC)/mbcTotalInflow;
            }

            scalar mbcErrorCeliaTimestep = -999;
            scalar mbcTimestepError = -999;
            if (mbcTimestepInflow != 0.0) {
                mbcErrorCeliaTimestep = (mbcOldTimestepTWC - mbcTimestepTWC)/mbcTimestepInflow;
                mbcTimestepError = mbcOldTimestepTWC - mbcTimestepTWC - mbcTimestepInflow;
            }

            scalar mbcTimestepErrorRelative = -999;
            if (mbcOldTimestepTWC != 0.0) {
                mbcTimestepErrorRelative = mbcTimestepError / mbcOldTimestepTWC;
            }
            
            mbcOldTimestepTWC = mbcTimestepTWC;

            if (Pstream::master())
            {
                LogMBC << runTime.value() << ";" << mbcTimestepError << ";" << mbcTimestepErrorRelative << ";" << mbcError << ";" << mbcErrorRelative << ";" << mbcErrorCelia << ";" << mbcErrorCeliaTimestep << ";" << mbcDiffTWC << ";" << mbcDiffRelativeTWC << ";" << mbcTotalInflow << ";" << mbcDiffTWCTimestep << ";" << mbcTimestepInflow<< ";" <<config.algorithm.getMbcRowInfo().c_str()<< nl;
                if (is_force_flush) {
                    LogMBC.flush();
                }
            }
        }
        // STOP - MBC timestep evaluation

#ifdef ADAPTIVE_MESH
        magGradH = mag(fvc::grad(h));
#endif

#include "write.H"

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << nl;

        // if (getFieldValueAtHeight(h, mesh, 0.025277778) < min(retentionModel->getShpHa()).value())
        // {   
        //     Info<<"Simulation stopped due to achieving h_a potential."<<nl;
        //     return 0;
        // } 
        
        runTime.printExecutionTime(Info);
    }

    //SIMPLE loop end
    Info << "Simulation Finished OK" << nl;

    return 0;
}
