#include "readGravitationalAcceleration.H"

using namespace Soil::RetentionModels;

// determine component aligned with non zero gravity vector
const bool hasX = mag(g.component(vector::X)).value() != 0.0;
const bool hasY = mag(g.component(vector::Y)).value() != 0.0;
const bool hasZ = mag(g.component(vector::Z)).value() != 0.0;

if (hasX + hasY + hasZ != 1)
{
    Info << "Gravity vector is not aligned with any of the coordinate axes. Exiting now." << endl;
    exit(0);
}

int gNonZeroComp = 0;
if (hasX)
    gNonZeroComp = vector::X;
if (hasY)
    gNonZeroComp = vector::Y;
if (hasZ)
    gNonZeroComp = vector::Z;

volScalarField h(
    IOobject(
        "h",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE),
    mesh);

volScalarField Theta(
    IOobject(
        "Theta",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("Theta", dimensionSet(0, 0, 0, 0, 0, 0, 0), scalar(0)));

volScalarField K(
    IOobject(
        "K",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("K", dimensionSet(0, 1, -1, 0, 0, 0, 0), 1));

volScalarField Cv(
    IOobject(
        "Cv",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("Cv", dimensionSet(0, -1, 0, 0, 0, 0, 0), 1));



#ifdef ADAPTIVE_MESH
volScalarField magGradH(
    IOobject(
        "magGradH",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    mag(fvc::grad(h)));
#endif

volScalarField region(
    IOobject(
        "region",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("region", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0));

volScalarField sources(
    IOobject(
        "sources",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("sources", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0));

word retentionModelString = word(transportProperties.lookup("swrc_model"));
retentionDataInterface *retentionModel = nullptr;

DEBUG("Retention model: " << retentionModelString);
if (config.processes.isWaterFilmFlow())
{

    if (retentionModelString == "vg")
    {
        if (config.processes.getWaterVaporFlow() == "vapor_std")
        {
            retentionModel = new retentionDataFilmVaporStdVg(transportProperties, region, mesh, runTime, config);
        } 
        else if (config.processes.getWaterVaporFlow() == "vapor_no_storage")
        {    
            retentionModel = new retentionDataFilmVaporNsVg(transportProperties, region, mesh, runTime, config);
        } 
        else 
        {
            retentionModel = new retentionDataFilmVg(transportProperties, region, mesh, runTime, config);
        }
    }

    if (retentionModelString == "ks")
    {
        retentionModel = new retentionDataFilmKs(transportProperties, region, mesh, runTime, config);
    }

    if (retentionModelString == "bc")
    {
        Info << "BC retention model can't be used with water film flow model. Exiting now." << endl;
        exit(0);
    }

    if (retentionModelString == "exp")
    {
        Info << "EXP retention model can't be used with water film flow model. Exiting now." << endl;
        exit(0);
    }

    if (retentionModelString == "hvk")
    {
        Info << "HVK retention model can't be used with water film flow model. Exiting now." << endl;
        exit(0);
    }
}
else
{
    if (retentionModelString == "vg")
    {
        retentionModel = new retentionDataVg(transportProperties, region, mesh, runTime, config);
    }

    if (retentionModelString == "ks")
    {
        retentionModel = new retentionDataKs(transportProperties, region, mesh, runTime, config);
    }

    if (retentionModelString == "bc")
    {
        retentionModel = new retentionDataBc(transportProperties, region, mesh, runTime, config);
    }

    if (retentionModelString == "exp")
    {
        retentionModel = new retentionDataExp(transportProperties, region, mesh, runTime, config);
    }

    if (retentionModelString == "hvk")
    {
        retentionModel = new retentionDataHvk(transportProperties, region, mesh, runTime, config);
    }
}


if (config.processes.getWaterVaporFlow() != "none" && !config.processes.isWaterFilmFlow())
{
    Info<<"Error! vapor_flow may be enabled only if film_flow is also enabled.";
    exit(0);
};

if (retentionModel == nullptr)
{
    Info << "Missing retention model declaration in dictionary \"transportProperties\". Exiting now." << endl;
    exit(0);
}

Theta = retentionModel->Theta(h, Theta);
K = retentionModel->Kh(h, K);
Cv = retentionModel->Cv(h, Cv);


if (config.algorithm.getShpDebug()) {
    debugShp(mesh, runTime, retentionModel);
    Info<<"Exiting after saving SHP info."<<nl;
    exit(0);
}


#ifdef VSSFOAM_DEBUG

volVectorField fluxWater(
    IOobject(
        "fluxWater",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedVector("fluxWater", dimensionSet(0, 1, -1, 0, 0, 0, 0), vector(0, 0, 0)));
fluxWater = K * (fvc::grad(h) + g / mag(g));

volScalarField fluxWaterX (
    IOobject(
        "fluxWaterX",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    fluxWater.component(vector::X));

volScalarField fluxWaterY (
    IOobject(
        "fluxWaterY",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    fluxWater.component(vector::Y));

volScalarField fluxWaterZ (
    IOobject(
        "fluxWaterZ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    fluxWater.component(vector::Z));


volScalarField hTotal(
    IOobject(
        "hTotal",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("a", dimensionSet(0, 1, 0, 0, 0, 0, 0), 0));
hTotal = h - mesh.C().component(gNonZeroComp);

volVectorField gradH(fvc::grad(h));

volScalarField gradHx(
    IOobject(
        "gradHx",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    gradH.component(vector::X));

volScalarField gradHy(
    IOobject(
        "gradHy",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    gradH.component(vector::Y));

volScalarField gradHz(
    IOobject(
        "gradHz",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    gradH.component(vector::Z));

volScalarField Kfilm(
    IOobject(
        "Kfilm",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("Kfilm", dimensionSet(0, 1, -1, 0, 0, 0, 0), 1));

#endif
