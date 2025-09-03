if (runTime.writeTime())
{
    runTime.write();

#ifdef VSSFOAM_DEBUG
    gradH = fvc::grad(h);

    gradHx = gradH.component(vector::X);
    gradHx.write();

    gradHy = gradH.component(vector::Y);
    gradHy.write();

    gradHz = gradH.component(vector::Z);
    gradHz.write();

    fluxWater = K * (fvc::grad(h) + g / mag(g));
    fluxWater.write();

    fluxWaterX = fluxWater.component(vector::X);
    fluxWaterX.write();

    fluxWaterY = fluxWater.component(vector::Y);
    fluxWaterY.write();

    fluxWaterZ = fluxWater.component(vector::Z);
    fluxWaterZ.write();

    hTotal = h - mesh.C().component(gNonZeroComp);
    hTotal.write();

    //Kfilm = retentionModel->Klf(h);
    //Kfilm.write();

    gradK.write();

    retentionModel->write();

#endif

}
