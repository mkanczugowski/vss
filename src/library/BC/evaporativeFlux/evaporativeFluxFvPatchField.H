#ifndef evaporativeFluxFvPatchField_H
#define evaporativeFluxFvPatchField_H

#include "fvPatchFields.H"
#include "fixedFluxFvPatchField.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    /*---------------------------------------------------------------------------*\
                     Class fixedFluxFvPatchField Declaration
    \*---------------------------------------------------------------------------*/
    class evaporativeFluxFvPatchField
        : public fixedFluxFvPatchField // TODO: Use zeroGradientFvPatchField as the parent class
    {
        // Private data
        scalar wind_speed_;             // Wind speed
        scalar surface_roughness_; // Surface roughness length
        scalar reference_level_;   // Reference level, the Z(vertical) coordinate of the reference level (to calculate the distance from soil surface)
        scalar air_temperature_;        // Air temperature at the reference height (in Celsius)
        scalar soil_temperature_drop_;  // Assumed soil surface temperature drop in relation to air temperature (in Celsius)
        bool is_first_run = true;

    public:
        //- Runtime type information
        TypeName("evaporativeFlux");

        // Constructors
        //- Construct from patch and internal field
        evaporativeFluxFvPatchField(
            const fvPatch &,
            const DimensionedField<scalar, volMesh> &);

        //- Construct from patch, internal field and dictionary
        evaporativeFluxFvPatchField(
            const fvPatch &,
            const DimensionedField<scalar, volMesh> &,
            const dictionary &);

        //- Construct by mapping the given evaporativeFluxFvPatchField
        //  onto a new patch
        evaporativeFluxFvPatchField(
            const evaporativeFluxFvPatchField &,
            const fvPatch &,
            const DimensionedField<scalar, volMesh> &,
            const fvPatchFieldMapper &);

        //- Construct as copy
        evaporativeFluxFvPatchField(
            const evaporativeFluxFvPatchField &);

        //- Construct as copy setting internal field reference
        evaporativeFluxFvPatchField(
            const evaporativeFluxFvPatchField &,
            const DimensionedField<scalar, volMesh> &);

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>(
                new evaporativeFluxFvPatchField(*this));
        };

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone(
            const DimensionedField<scalar, volMesh> &iF) const
        {
            return tmp<fvPatchScalarField>(
                new evaporativeFluxFvPatchField(*this, iF));
        };

        void rmap(
            const fvPatchScalarField &ptf,
            const labelList &addr);

        tmp<scalarField> calculateGradient(void) const;

        virtual void write(Ostream &) const;
    };
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
