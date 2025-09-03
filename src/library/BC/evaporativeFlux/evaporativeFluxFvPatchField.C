#include "evaporativeFluxFvPatchField.H"
#include "dictionary.H"
#include "fvCFD.H"
#include "gravityMeshObject.H"
#include "addToRunTimeSelectionTable.H"

#undef VSSFOAM_DEBUG

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
    //- Construct from patch and internal field
    evaporativeFluxFvPatchField::evaporativeFluxFvPatchField(
        const fvPatch &p,
        const DimensionedField<scalar, volMesh> &iF)
        : fixedFluxFvPatchField(p, iF),
          wind_speed_(0),
          surface_roughness_(0),
          reference_level_(0),
          air_temperature_(0),
          soil_temperature_drop_(0)
    {
    }

    //- Construct from patch, internal field and dictionary
    evaporativeFluxFvPatchField::evaporativeFluxFvPatchField(
        const fvPatch &p,
        const DimensionedField<scalar, volMesh> &iF,
        const dictionary &dict)
        : fixedFluxFvPatchField(p, iF, dict),
          wind_speed_(dict.get<scalar>("wind_speed")),
          surface_roughness_(dict.get<scalar>("surface_roughness")),
          reference_level_(dict.get<scalar>("reference_level")),
          air_temperature_(dict.get<scalar>("air_temperature")),
          soil_temperature_drop_(dict.get<scalar>("soil_temperature_drop"))
    {
        Info << "evaporativeFluxFvPatchField" << nl;
    }

    //- Construct by mapping the given evaporativeFluxFvPatchField
    //  onto a new patch
    evaporativeFluxFvPatchField::evaporativeFluxFvPatchField(
        const evaporativeFluxFvPatchField &ptf,
        const fvPatch &p,
        const DimensionedField<scalar, volMesh> &iF,
        const fvPatchFieldMapper &mapper)
        : fixedFluxFvPatchField(ptf, p, iF, mapper),
          wind_speed_(ptf.wind_speed_),
          surface_roughness_(ptf.surface_roughness_),
          reference_level_(ptf.reference_level_),
          air_temperature_(ptf.air_temperature_),
          soil_temperature_drop_(ptf.soil_temperature_drop_)
    {
    }

    //- Construct as copy
    evaporativeFluxFvPatchField::evaporativeFluxFvPatchField(
        const evaporativeFluxFvPatchField &ptf)
        : fixedFluxFvPatchField(ptf),
          wind_speed_(ptf.wind_speed_),
          surface_roughness_(ptf.surface_roughness_),
          reference_level_(ptf.reference_level_),
          air_temperature_(ptf.air_temperature_),
          soil_temperature_drop_(ptf.soil_temperature_drop_)
    {
    }

    //- Construct as copy setting internal field reference
    evaporativeFluxFvPatchField::evaporativeFluxFvPatchField(
        const evaporativeFluxFvPatchField &ptf,
        const DimensionedField<scalar, volMesh> &iF)
        : fixedFluxFvPatchField(ptf, iF),
          wind_speed_(ptf.wind_speed_),
          surface_roughness_(ptf.surface_roughness_),
          reference_level_(ptf.reference_level_),
          air_temperature_(ptf.air_temperature_),
          soil_temperature_drop_(ptf.soil_temperature_drop_)
    {
    }

    // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
    //- Reverse map the given fvPatchField onto this fvPatchField
    void evaporativeFluxFvPatchField::rmap(
        const fvPatchScalarField &ptf,
        const labelList &addr)
    {
        // fvPatchScalarField::rmap(ptf, addr);

        // const evaporativeFluxFvPatchField &fgptf =
        //     refCast<const evaporativeFluxFvPatchField>(ptf);

        // flux().rmap(fgptf.flux(), addr);

        // this->initiateFluxNormVert();
    }

    tmp<scalarField> evaporativeFluxFvPatchField::calculateGradient(void) const
    {
        const word patchName = patch().name();
        tmp<scalarField> K = patch().lookupPatchField<volScalarField>("K");
        tmp<scalarField> h = patch().lookupPatchField<volScalarField>("h");
        tmp<scalarField> Theta = patch().lookupPatchField<volScalarField>("Theta");
        tmp<scalarField> shp_ret_th_s = patch().lookupPatchField<volScalarField>("shp_ret_th_s");
        const uniformDimensionedVectorField &g_vector = meshObjects::gravity::New(this->db().time());
        int vertical_component = -1;
        if (g_vector.component(0).value() != 0.0)
            vertical_component = 0;
        if (g_vector.component(1).value() != 0.0)
            vertical_component = 1;
        if (g_vector.component(2).value() != 0.0)
            vertical_component = 2;
        const scalar rho_lw = 1000;
        const scalar g = Foam::mag(g_vector.value());
        vectorField patch_cf = patch().Cf(); //.component(vertical_component);
        scalarField z_patch = patch_cf.component(vertical_component);
        scalarField z_r = reference_level_ - z_patch;
        const scalar T_ar = air_temperature_ + 273.15;
        const scalar T_srf = air_temperature_ + 273.15 - soil_temperature_drop_;
        const scalar rho_sv_T_srf = Foam::exp(31.3716 - 6014.79 / T_srf - 7.92495e-3 * T_srf) / (rho_lw * T_srf);
        const scalar rho_sv_T_ar = Foam::exp(31.3716 - 6014.79 / T_ar - 7.92495e-3 * T_ar) / (rho_lw * T_ar);
        const scalar Rsv = 461.52;
        scalarField rho_vs = rho_sv_T_srf * Foam::exp((h * g) / (Rsv * T_srf)); // Saito 2006
        const scalar rho_va = rho_sv_T_ar;
        tmp<scalarField> R_i = g * (z_r - surface_roughness_) * (soil_temperature_drop_) / (T_ar * pow(wind_speed_, 2));
        tmp<scalarField> S_t = 1 / (1 - 10 * R_i);
        scalarField r_v = Foam::pow(Foam::log(z_r / surface_roughness_), 2) * S_t / (0.16 * wind_speed_);
        //scalarField r_s = -805+4140*(ThetaSat-Theta); //Saito 2006 - not working correctly
        //scalarField r_s = 10 * Foam::max(One, Foam::exp(35.36 * (0.15 - Theta))); 
        scalarField r_s = 33.5+3.5*Foam::pow(shp_ret_th_s/Theta, 2.3); //Shu Fen Sun 1982 - gives reasonable results
        scalarField E = (rho_vs - rho_va) / (r_v + r_s);
        scalarField q_le = E / rho_lw;
        scalarField tmp = q_le / K;
#ifdef VSSFOAM_DEBUG
        const scalar t = this->db().time().timeOutputValue();
        Info << "EVAPINFO;" << t << ";" << T_ar << ";" << T_srf << ";" << rho_sv_T_srf << ";" << rho_sv_T_ar << ";" << max(rho_vs) << ";" << rho_va << ";" << max(r_v) << ";" << max(r_s) << ";" << max(shp_ret_th_s) << ";" << max(Theta) << ";" << max(E) << ";" << max(q_le) << ";" << max(K) << nl;
#endif
        return (gradient() + tmp);
    }

    void evaporativeFluxFvPatchField::write(Ostream &os) const
    {
        zeroFluxFvPatchField::write(os);
        flux().writeEntry("flux", os);
        os.writeEntry("wind_speed", wind_speed_);
        os.writeEntry("surface_roughness", surface_roughness_);
        os.writeEntry("reference_level", reference_level_);
        os.writeEntry("air_temperature", air_temperature_);
        os.writeEntry("soil_temperature_drop", soil_temperature_drop_);
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam

namespace Foam
{
    makePatchTypeField(
        fvPatchScalarField,
        evaporativeFluxFvPatchField);
}; // End namespace Foam

// ************************************************************************* //
