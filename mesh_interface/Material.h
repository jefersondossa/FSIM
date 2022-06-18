#ifndef MATERIAL_H
#define MATERIAL_H

enum class MaterialType
{
	ELASTIC_SOLID,
	ELASTIC_INCOMPRESSIBLE_SOLID,
	NEWTONIAN_INCOMPRESSIBLE_FLUID
};

enum PlaneAnalysis
{
	PLANE_STRESS,
	PLANE_STRAIN
};

class Material
{
public:
	Material(const double& density, const MaterialType& type, const PlaneAnalysis& planeAnalysis = PLANE_STRESS);

	virtual ~Material() = 0;

	double getDensity() const;

	MaterialType getType() const;

	PlaneAnalysis getPlaneAnalysis() const;

	void setDensity(const double& density);

	void setType(const MaterialType& type);

	void setPlaneAnalysis(const PlaneAnalysis& planeAnalysis);

	virtual void getPlaneStressTensor(const double E[3], const double dx_dy[2][2], double S[3]) = 0;

	virtual void getPlaneStressTensorDerivative(const double dE_dy[3], const double dx_dy[2][2], double dS_dy[3]) = 0;

	virtual void getStressTensor(const double E[6], const double dx_dy[3][3], double S[6]) = 0;

	virtual void getStressTensorDerivative(const double dE_dy[6], const double dx_dy[3][3], double dS_dy[6]) = 0;

protected:
	double density_;
	MaterialType type_;
	PlaneAnalysis planeAnalysis_;
};

class ElasticSolid : public Material
{
	public:
		ElasticSolid(const double& young, const double& poisson, const double& density = 0.0);

		~ElasticSolid();

		double getYoung() const;

		double getPoisson() const;

		void setYoung(const double& young);

		void setPoisson(const double& poisson);

		void getPlaneStressTensor(const double E[3], const double dx_dy[2][2], double S[3]) override;

		void getPlaneStressTensorDerivative(const double dE_dy[3], const double dx_dy[2][2], double dS_dy[3]) override;

		void getStressTensor(const double E[6], const double dx_dy[3][3], double S[6]) override;

		void getStressTensorDerivative(const double dE_dy[6], const double dx_dy[3][3], double dS_dy[6]) override;

	private:
		double young_;
		double poisson_;
};

class NewtonianFluid : public Material
{
	public:
		NewtonianFluid(const double& viscosity, const double& density);

		~NewtonianFluid();

		double getViscosity() const;

		void setViscosity(const double& viscosity);

		virtual void getPlaneStressTensor(const double dE_dt[3], const double dx_dy[2][2], double S[3]) override;

		virtual void getPlaneStressTensorDerivative(const double dE_dtdy[3], const double dx_dy[2][2], double dS_dy[3]) override;

		virtual void getStressTensor(const double dE_dt[6], const double dx_dy[3][3], double S[6]) override;

		virtual void getStressTensorDerivative(const double dE_dtdy[6], const double dx_dy[3][3], double dS_dy[6]) override;

	private:
		double viscosity_;
};

#endif



















