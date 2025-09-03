#include "retentionData.h"
#include "soilPhysics.h"

using namespace Soil::RetentionModels;
using namespace Soil::Physics;

namespace Soil::RetentionModels
{

	retentionData::retentionData(IOdictionary &transportProperties, volScalarField &region, fvMesh &mesh, Time &runTime, simConfig &config)
		: retentionDataInterface(mesh, runTime, config),
		  soil_temp(IOobject("soil_temperature", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("soil_temperature", dimensionSet(0, 0, 0, 1, 0, 0, 0), 293.15)),
		  shp_cv_res(IOobject("shp_cv_res", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("shp_cv_res", dimensionSet(0, -1, 0, 0, 0, 0, 0), 0.0)),
		  shp_mual_se_l(IOobject("shp_mual_se_l", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("shp_mual_se_l", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.5)),
		  shp_dne_tau(IOobject("shp_dne_tau", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("shp_dne_tau", dimensionSet(0, 0, 1, 0, 0, 0, 0), 1)),
		  shp_mual_k_sat(IOobject("shp_mual_k_sat", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("shp_mual_k_sat", dimensionSet(0, 1, -1, 0, 0, 0, 0), 1)),
		  shp_ret_th_s(IOobject("shp_ret_th_s", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("shp_ret_th_s", dimensionSet(0, 0, 0, 0, 0, 0, 0), 1)),
		  shp_ret_th_r(IOobject("shp_ret_th_r", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("shp_ret_th_r", dimensionSet(0, 0, 0, 0, 0, 0, 0), 1)),
		  shp_ret_total_por(IOobject("shp_ret_total_por", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("shp_ret_total_por", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0))
	{
		char buffer[1000];
		no_of_regions = static_cast<int>(transportProperties.lookupOrDefault<scalar>("no_of_regions", 1.0));
		Info << "nr. of regions:" << no_of_regions << nl;

		
		scalar soil_temperature = static_cast<scalar>(transportProperties.lookupOrDefault<scalar>("soil_temperature", 293.15));

		word dneModel = config.processes.getDneModel();

		if (no_of_regions > 0)
		{
			shp_mual_se_l_vector.reserve(no_of_regions);
			for (int i = 0; i < no_of_regions; i++)
			{
				sprintf(buffer, "shp_mual_se_l_%d", i);
				std::string key_name(buffer);
				shp_mual_se_l_vector.push_back(dimensionedScalar(key_name.c_str(), dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.5, transportProperties));
			}

			if (dneModel == "ross_smettem")
			{
				Info << "init shp_dne_tau" << nl;
				shp_dne_tau_vector.reserve(no_of_regions);
				for (int i = 0; i < no_of_regions; i++)
				{
					sprintf(buffer, "shp_dne_tau_%d", i);
					std::string key_name(buffer);
					shp_dne_tau_vector.push_back(dimensionedScalar(key_name.c_str(), dimensionSet(0, 0, 1, 0, 0, 0, 0), 1.0, transportProperties));
					Info << dimensionedScalar(key_name.c_str(), transportProperties) << nl;
				}
			}

			shp_mual_k_sat_vector.reserve(no_of_regions);
			for (int i = 0; i < no_of_regions; i++)
			{
				sprintf(buffer, "shp_mual_k_sat_%d", i);
				std::string key_name(buffer);
				shp_mual_k_sat_vector.push_back(dimensionedScalar(key_name.c_str(), dimensionSet(0, 1, -1, 0, 0, 0, 0), transportProperties));
			}

			shp_ret_th_s_vector.reserve(no_of_regions);
			for (int i = 0; i < no_of_regions; i++)
			{
				sprintf(buffer, "shp_ret_th_s_%d", i);
				std::string key_name(buffer);
				shp_ret_th_s_vector.push_back(dimensionedScalar(key_name.c_str(), dimensionSet(0, 0, 0, 0, 0, 0, 0), transportProperties));
			}

			shp_ret_th_r_vector.reserve(no_of_regions);
			for (int i = 0; i < no_of_regions; i++)
			{
				sprintf(buffer, "shp_ret_th_r_%d", i);
				std::string key_name(buffer);
				shp_ret_th_r_vector.push_back(dimensionedScalar(key_name.c_str(), dimensionSet(0, 0, 0, 0, 0, 0, 0), transportProperties));
			}

			shp_ret_total_por_vector.reserve(no_of_regions);
			for (int i = 0; i < no_of_regions; i++)
			{
				sprintf(buffer, "shp_ret_total_por_%d", i);
				std::string key_name(buffer);
				shp_ret_total_por_vector.push_back(dimensionedScalar(key_name.c_str(), dimensionSet(0, 0, 0, 0, 0, 0, 0), shp_ret_th_s_vector[i].value(), transportProperties));
			}

			forAll(region, index)
			{
				if (dneModel == "ross_smettem")
				{
					shp_dne_tau[index] = shp_dne_tau_vector[static_cast<int>(region[index])].value();
					Info << shp_dne_tau[index] << nl;
				};
				shp_mual_se_l[index] = shp_mual_se_l_vector[static_cast<int>(region[index])].value();
				shp_mual_k_sat[index] = shp_mual_k_sat_vector[static_cast<int>(region[index])].value();
				shp_ret_th_s[index] = shp_ret_th_s_vector[static_cast<int>(region[index])].value();
				shp_ret_th_r[index] = shp_ret_th_r_vector[static_cast<int>(region[index])].value();
				shp_ret_total_por[index] = shp_ret_total_por_vector[static_cast<int>(region[index])].value();
				shp_cv_res[index] = config.algorithm.getShpCvRes();
			}

			forAll(region.boundaryFieldRef(), ipatch)
			{
				forAll(region.boundaryFieldRef()[ipatch], index)
				{
					if (dneModel == "ross_smettem")
					{
						shp_dne_tau.boundaryFieldRef()[ipatch][index] = shp_dne_tau_vector[static_cast<int>(region.boundaryField()[ipatch][index])].value();
					};
					shp_mual_se_l.boundaryFieldRef()[ipatch][index] = shp_mual_se_l_vector[static_cast<int>(region.boundaryField()[ipatch][index])].value();
					shp_mual_k_sat.boundaryFieldRef()[ipatch][index] = shp_mual_k_sat_vector[static_cast<int>(region.boundaryField()[ipatch][index])].value();
					shp_ret_th_s.boundaryFieldRef()[ipatch][index] = shp_ret_th_s_vector[static_cast<int>(region.boundaryField()[ipatch][index])].value();
					shp_ret_th_r.boundaryFieldRef()[ipatch][index] = shp_ret_th_r_vector[static_cast<int>(region.boundaryField()[ipatch][index])].value();
					shp_ret_total_por.boundaryFieldRef()[ipatch][index] = shp_ret_total_por_vector[static_cast<int>(region.boundaryField()[ipatch][index])].value();
					shp_cv_res.boundaryFieldRef()[ipatch][index] = config.algorithm.getShpCvRes();
					soil_temp.boundaryFieldRef()[ipatch][index] = soil_temperature;
				}
			}
		}
		else
		{
			if (dneModel == "ross_smettem")
			{
				shp_dne_tau = volScalarField(IOobject("shp_dne_tau", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);
			};
			shp_mual_se_l = volScalarField(IOobject("shp_mual_se_l", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);
			shp_mual_k_sat = volScalarField(IOobject("shp_mual_k_sat", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);
			shp_ret_th_s = volScalarField(IOobject("shp_ret_th_s", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);
			shp_ret_th_r = volScalarField(IOobject("shp_ret_th_r", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);
			shp_ret_total_por = volScalarField(IOobject("shp_ret_total_por", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);
		}
	}

	void retentionData::errInfo(void)
	{
		Info << "Missing retention model declaration in dictionary \"transportProperties\". Exiting now." << endl;
		exit(0);
	}
} // namespace Soil::RetentionModelsF
