#define _USE_MATH_DEFINES
#include <cmath>

#include "GreenCoordinates.h"

#include <iostream>
#include <vector>
#include <Eigen/Geometry>

/// Taken from https://github.com/superboubek/QMVC

template< class point_t >
point_t triMeanVector(point_t const& p, point_t const& p0, point_t const& p1, point_t const& p2)
{
	point_t pttri[3] = { p0,p1,p2 };
	point_t utri[3];
	for (unsigned int i = 0; i < 3; ++i) {
		utri[i] = (pttri[i] - p);
		if (utri[i].stableNorm() < 0.0000000000000001)
		{
			utri[i] = { 0, 0, 0 };
		}
		else
		{
			utri[i] = utri[i].normalized();
		}
		assert(utri[i] == utri[i]);
	}

	double thetatri[3];
	for (unsigned int i = 0; i < 3; ++i) {
		assert(!std::isnan((utri[(i) % 3] - utri[(i + 1) % 3]).norm()));
		thetatri[i] = 2.0 * std::asin(std::min<double>(1.0, std::max<double>(-1.0, (utri[(i) % 3] - utri[(i + 1) % 3]).norm() / 2.0)));
		assert(!std::isnan(thetatri[i]));
	}

	point_t Ntri[3];
	for (unsigned int i = 0; i < 3; ++i) {
		Ntri[i] = (pttri[(i + 1) % 3] - p).cross(pttri[(i) % 3] - p);
		assert(Ntri[i] == Ntri[i]);
	}

	auto Ntri_norm_0 = Ntri[0].stableNorm();
	if (Ntri_norm_0 < 0.0000000000000001)
	{
		Ntri[0] = { 0, 0, 0 };
	}
	else
	{
		Ntri[0] = Ntri[0].normalized();
	}
	auto Ntri_norm_1 = Ntri[1].stableNorm();
	if (Ntri_norm_1 < 0.0000000000000001)
	{
		Ntri[1] = { 0, 0, 0 };
	}
	else
	{
		Ntri[1] = Ntri[1].normalized();
	}
	auto Ntri_norm_2 = Ntri[2].stableNorm();
	if (Ntri_norm_2 < 0.0000000000000001)
	{
		Ntri[2] = { 0, 0, 0 };
	}
	else
	{
		Ntri[2] = Ntri[2].normalized();
	}
	point_t m_tri = -0.5 * (thetatri[0] * Ntri[0] + thetatri[1] * Ntri[1] + thetatri[2] * Ntri[2]);
	return m_tri;
}



template< class point_t >
point_t quadMeanVector(point_t const& p, point_t const& p0, point_t const& p1, point_t const& p2, point_t const& p3)
{
	point_t mv(0, 0, 0);
	point_t t0, t1, t2, pc = (p0 + p1 + p2 + p3) / 4.0;
	t0 = p0; t1 = p1; t2 = pc;
	{
		mv += triMeanVector(p, t0, t1, t2);
	}
	t0 = p1; t1 = p2; t2 = pc;
	{
		mv += triMeanVector(p, t0, t1, t2);
	}
	t0 = p2; t1 = p3; t2 = pc;
	{
		mv += triMeanVector(p, t0, t1, t2);
	}
	t0 = p3; t1 = p0; t2 = pc;
	{
		mv += triMeanVector(p, t0, t1, t2);
	}
	return mv;
}

template< class point_t >
double computeClosestPointOnLineFromLine(point_t const& A, point_t const& B, point_t const& eta, point_t const& d) {
	point_t const& e = B - A;
	double eDotD = e.dot(d);
	double dDotD = d.dot(d);
	double eDotE = e.dot(e);

	double alpha = eDotE - eDotD * eDotD / dDotD;
	double beta = (A - eta).dot(e) - (A - eta).dot(d) * eDotD / dDotD;

	return -beta / alpha;
}

template< class point_t >
void computeClosestPointOnQuadFromLine(point_t const& q0, point_t const& q1, point_t const& q2, point_t const& q3,
	point_t const& origin, point_t const& direction,
	double& uProj, double& vProj, bool ForceIntoUnitSquare = true) {
	double epsilon = 0.0;

	double aV = ((q3 - q0).cross(q2 - q1)).dot(direction);
	double aU = ((q1 - q0).cross(q2 - q3)).dot(direction);

	if (std::fabs(aV) > std::fabs(aU)) {
		// solve for v :
		double a = ((q3 - q0).cross(q2 - q1)).dot(direction);
		double b = ((q0 - origin).cross(q2 - q1) + (q3 - q0).cross(q1 - origin)).dot(direction);
		double c = ((q0 - origin).cross(q1 - origin)).dot(direction);
		if (std::fabs(a) > epsilon) {
			double Delta = b * b - 4 * a * c;
			if (Delta >= 0.0) {
				double vPlus = (-b + std::sqrt(Delta)) / (2 * a);
				double vMinus = (-b - std::sqrt(Delta)) / (2 * a);
				if (std::fabs(vMinus - 0.5) < std::fabs(vPlus - 0.5))
					vProj = vMinus;
				else
					vProj = vPlus;
			}
		}
		else {
			vProj = -c / b;
		}

		if (ForceIntoUnitSquare) vProj = std::min<double>(1.0, std::max<double>(0.0, vProj));

		point_t AV = (1.0 - vProj) * q0 + vProj * q3;
		point_t BV = (1.0 - vProj) * q1 + vProj * q2;

		uProj = computeClosestPointOnLineFromLine(AV, BV, origin, direction);
		if (ForceIntoUnitSquare) uProj = std::min<double>(1.0, std::max<double>(0.0, uProj));
	}
	else {
		// solve for u :
		double a = ((q1 - q0).cross(q2 - q3)).dot(direction);
		double b = ((q0 - origin).cross(q2 - q3) + (q1 - q0).cross(q3 - origin)).dot(direction);
		double c = ((q0 - origin).cross(q3 - origin)).dot(direction);
		if (std::fabs(a) > epsilon) {
			double Delta = b * b - 4 * a * c;
			if (Delta >= 0.0) {
				double uPlus = (-b + std::sqrt(Delta)) / (2 * a);
				double uMinus = (-b - std::sqrt(Delta)) / (2 * a);
				if (std::fabs(uMinus - 0.5) < std::fabs(uPlus - 0.5))
					uProj = uMinus;
				else
					uProj = uPlus;
			}
		}
		else {
			uProj = -c / b;
		}

		if (ForceIntoUnitSquare) uProj = std::min<double>(1.0, std::max<double>(0.0, uProj));

		point_t AU = (1.0 - uProj) * q0 + uProj * q1;
		point_t BU = (1.0 - uProj) * q3 + uProj * q2;

		vProj = computeClosestPointOnLineFromLine(AU, BU, origin, direction);
		if (ForceIntoUnitSquare) vProj = std::min<double>(1.0, std::max<double>(0.0, vProj));
	}
}

template< class float_t, class point_t >
bool computeUnnormalizedMVCForOneTriangle(
	point_t const& eta,
	point_t* tri_vertices,
	float_t* w_weights)
{
	typedef double   T;

	T epsilon = 0.000000001;

	T d[3];
	point_t u[3];

	for (unsigned int v = 0; v < 3; ++v)
	{
		d[v] = (eta - tri_vertices[v]).norm();
		u[v] = (tri_vertices[v] - eta) / d[v];
	}

	T l[3]; T theta[3];

	{
		// the Norm is CCW :
		for (unsigned int i = 0; i <= 2; ++i) {
			l[i] = (u[(i + 1) % 3] - u[(i + 2) % 3]).norm();
		}

		for (unsigned int i = 0; i <= 2; ++i) {
			theta[i] = 2.0 * asin(l[i] / 2.0);
		}

		T determinant = (tri_vertices[0] - eta).dot(
			(tri_vertices[1] - tri_vertices[0]).cross(tri_vertices[2] - tri_vertices[0]));
		T sqrdist = determinant * determinant / (4 * (tri_vertices[1] - tri_vertices[0]).cross(
			tri_vertices[2] - tri_vertices[0]).squaredNorm());
		T dist = std::sqrt(sqrdist);

		if (dist < epsilon)
		{
			// then the point eta lies on the support plane of the triangle
			T h = (theta[0] + theta[1] + theta[2]) / 2.0;
			if (M_PI - h < epsilon)
			{
				// then eta lies inside the triangle t
				return true;
			}
			else {
				w_weights[0] = w_weights[1] = w_weights[2] = 0.0;
				return false;
			}
		}

		point_t N[3];

		for (unsigned int i = 0; i < 3; ++i)
			N[i] = (tri_vertices[(i + 1) % 3] - eta).cross(tri_vertices[(i + 2) % 3] - eta);

		for (unsigned int i = 0; i <= 2; ++i)
		{
			w_weights[i] = 0.0;
			for (unsigned int j = 0; j <= 2; ++j)
				w_weights[i] += theta[j] * (N[i].dot(N[j])) / (2.0 * N[j].norm());

			w_weights[i] /= determinant;
		}
	}
	return false;
}

template< class point_t >
point_t smoothProjectInsideTet(point_t const& eta, point_t const& q0, point_t const& q1, point_t const& q2, point_t const& q3) {
	point_t proj(0, 0, 0);
	point_t tri[3];
	double sumWeights = 0.0;
	double phi[3];

	tri[0] = q0; tri[1] = q1; tri[2] = q2;
	{
		computeUnnormalizedMVCForOneTriangle(eta, tri, phi);
		sumWeights += fabs(phi[0]) + fabs(phi[1]) + fabs(phi[2]);
		proj += fabs(phi[0]) * tri[0] + fabs(phi[1]) * tri[1] + fabs(phi[2]) * tri[2];
	}

	tri[0] = q0; tri[1] = q1; tri[2] = q3;
	{
		computeUnnormalizedMVCForOneTriangle(eta, tri, phi);
		sumWeights += fabs(phi[0]) + fabs(phi[1]) + fabs(phi[2]);
		proj += fabs(phi[0]) * tri[0] + fabs(phi[1]) * tri[1] + fabs(phi[2]) * tri[2];
	}

	tri[0] = q0; tri[1] = q3; tri[2] = q2;
	{
		computeUnnormalizedMVCForOneTriangle(eta, tri, phi);
		sumWeights += fabs(phi[0]) + fabs(phi[1]) + fabs(phi[2]);
		proj += fabs(phi[0]) * tri[0] + fabs(phi[1]) * tri[1] + fabs(phi[2]) * tri[2];
	}

	tri[0] = q3; tri[1] = q1; tri[2] = q2;
	{
		computeUnnormalizedMVCForOneTriangle(eta, tri, phi);
		sumWeights += fabs(phi[0]) + fabs(phi[1]) + fabs(phi[2]);
		proj += fabs(phi[0]) * tri[0] + fabs(phi[1]) * tri[1] + fabs(phi[2]) * tri[2];
	}

	return proj / sumWeights;
}

template< class point_t >
point_t smoothProjectInsideTet_second(point_t const& eta, point_t const& q0, point_t const& q1, point_t const& q2, point_t const& q3) {
	point_t proj(0, 0, 0);
	point_t tri[3];
	double sumWeights = 0.0;
	double phi[3];

	tri[0] = q0; tri[1] = q1; tri[2] = q2;
	{
		computeUnnormalizedMVCForOneTriangle(eta, tri, phi);
		double dist_to_t = (tri[0] - eta).dot(((tri[1] - tri[0]).cross(tri[2] - tri[0])).normalized());
		if (fabs(dist_to_t) != 0.0) {
			phi[0] /= dist_to_t;
			phi[1] /= dist_to_t;
			phi[2] /= dist_to_t;
		}
		sumWeights += phi[0] + phi[1] + phi[2];
		proj += phi[0] * tri[0] + phi[1] * tri[1] + phi[2] * tri[2];
	}

	tri[0] = q0; tri[1] = q1; tri[2] = q3;
	{
		computeUnnormalizedMVCForOneTriangle(eta, tri, phi);
		double dist_to_t = (tri[0] - eta).dot(((tri[1] - tri[0]).cross(tri[2] - tri[0])).normalized());
		if (fabs(dist_to_t) != 0.0) {
			phi[0] /= dist_to_t;
			phi[1] /= dist_to_t;
			phi[2] /= dist_to_t;
		}
		sumWeights += phi[0] + phi[1] + phi[2];
		proj += phi[0] * tri[0] + phi[1] * tri[1] + phi[2] * tri[2];
	}

	tri[0] = q0; tri[1] = q3; tri[2] = q2;
	{
		computeUnnormalizedMVCForOneTriangle(eta, tri, phi);
		double dist_to_t = (tri[0] - eta).dot(((tri[1] - tri[0]).cross(tri[2] - tri[0])).normalized());
		if (fabs(dist_to_t) != 0.0) {
			phi[0] /= dist_to_t;
			phi[1] /= dist_to_t;
			phi[2] /= dist_to_t;
		}
		sumWeights += phi[0] + phi[1] + phi[2];
		proj += phi[0] * tri[0] + phi[1] * tri[1] + phi[2] * tri[2];
	}

	tri[0] = q3; tri[1] = q1; tri[2] = q2;
	{
		computeUnnormalizedMVCForOneTriangle(eta, tri, phi);
		double dist_to_t = (tri[0] - eta).dot(((tri[1] - tri[0]).cross(tri[2] - tri[0])).normalized());
		if (fabs(dist_to_t) != 0.0) {
			phi[0] /= dist_to_t;
			phi[1] /= dist_to_t;
			phi[2] /= dist_to_t;
		}
		sumWeights += phi[0] + phi[1] + phi[2];
		proj += phi[0] * tri[0] + phi[1] * tri[1] + phi[2] * tri[2];
	}

	return proj / sumWeights;
}

template<class point_t>
bool isInConvexHull(point_t const& p, point_t const& p0, point_t const& p1, point_t const& p2, point_t const& p3) {
	double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
	point_t t0, t1, t2;
	t0 = p0; t1 = p1; t2 = p3;
	{
		point_t nTri = (t1 - t0).cross(t2 - t0);
		s0 = (t0 - p).dot(nTri);
	}
	t0 = p1; t1 = p2; t2 = p3;
	{
		point_t nTri = (t1 - t0).cross(t2 - t0);
		s1 = (t0 - p).dot(nTri);
	}
	t0 = p0; t1 = p2; t2 = p1;
	{
		point_t nTri = (t1 - t0).cross(t2 - t0);
		s2 = (t0 - p).dot(nTri);
	}
	t0 = p0; t1 = p3; t2 = p2;
	{
		point_t nTri = (t1 - t0).cross(t2 - t0);
		s3 = (t0 - p).dot(nTri);
	}
	return (s0 >= 0.0 && s1 >= 0.0 && s2 >= 0.0 && s3 >= 0.0) || (s0 <= 0.0 && s1 <= 0.0 && s2 <= 0.0 && s3 <= 0.0);
}

template< class point_t >
void smoothProjectOnQuad(point_t const& eta, point_t* quad_vertices, double& uProject, double& vProject) {
	point_t pToProject = eta, m;
	bool pointIsStrictlyInConvexHull = true;

	if (!isInConvexHull(pToProject, quad_vertices[0], quad_vertices[1], quad_vertices[2], quad_vertices[3])) {
		{
			pToProject = smoothProjectInsideTet_second(eta, quad_vertices[0], quad_vertices[1], quad_vertices[2], quad_vertices[3]); // either smoothProjectInsideTet or smoothProjectInsideTet_second... try them
		}
	}

	{
		m = quadMeanVector(pToProject, quad_vertices[0], quad_vertices[1], quad_vertices[2], quad_vertices[3]).normalized();

		if (pointIsStrictlyInConvexHull)
			assert(m == m && "m.isnan() : bilMeanVector seems to be problematic for a point strictly inside the convex hull");
		else
			assert(m == m && "m.isnan() : bilMeanVector seems to be problematic for a point ON the convex hull");

		computeClosestPointOnQuadFromLine(quad_vertices[0], quad_vertices[1], quad_vertices[2], quad_vertices[3], pToProject, m, uProject, vProject);
	}
}

template< class point_t >
double get_signed_solid_angle(point_t const& a, point_t const& b, point_t const& c) {
	typedef double    T;
	T det = a.dot(b.cross(c));
	if (fabs(det) < 0.0000000001) // then you're on the limit case where you cover half the sphere
		return 2.0 * M_PI; // that is particularly shitty, because the sign is difficult to estimate...

	T al = a.norm(), bl = b.norm(), cl = c.norm();

	T div = al * bl * cl + a.dot(b) * cl + a.dot(c) * bl + b.dot(c) * al;
	T at = std::atan2(std::abs(det), div);
	if (at < 0) at += M_PI; // If det>0 && div<0 atan2 returns < 0, so add pi.
	T omega = 2.0 * at;

	if (det > 0.0) return omega;
	return -omega;
}

template< class float_t, class point_t >
void computePhiAndPsiForOneTriangle(point_t const& eta,
	point_t* tri_vertices, // an array of 3 points
	float_t* phi, // an array of 3 floats
	float_t& psi) {
	typedef double    T;
	point_t Nt = (tri_vertices[1] - tri_vertices[0]).cross(tri_vertices[2] - tri_vertices[0]);
	T NtNorm = Nt.norm();
	T At = NtNorm / 2.0;
	Nt /= NtNorm;

	psi = 0.0;
	for (unsigned int v = 0; v < 3; ++v) phi[v] = 0.0;

	point_t e[3];    T e_norm[3];   point_t e_normalized[3];    T R[3];    point_t d[3];    T d_norm[3];     T C[3];     point_t J[3];
	for (unsigned int v = 0; v < 3; ++v) e[v] = tri_vertices[v] - eta;
	for (unsigned int v = 0; v < 3; ++v) e_norm[v] = e[v].norm();
	for (unsigned int v = 0; v < 3; ++v) e_normalized[v] = e[v] / e_norm[v];

	T signed_solid_angle = get_signed_solid_angle(e_normalized[0], e_normalized[1], e_normalized[2]) / (4.f * M_PI);
	T signed_volume = (e[0].cross(e[1])).dot(e[2]) / 6.0;

	for (unsigned int v = 0; v < 3; ++v) R[v] = e_norm[(v + 1) % 3] + e_norm[(v + 2) % 3];
	for (unsigned int v = 0; v < 3; ++v) d[v] = tri_vertices[(v + 1) % 3] - tri_vertices[(v + 2) % 3];
	for (unsigned int v = 0; v < 3; ++v) d_norm[v] = d[v].norm();
	for (unsigned int v = 0; v < 3; ++v) C[v] = std::log((R[v] + d_norm[v]) / (R[v] - d_norm[v])) / (4.0 * M_PI * d_norm[v]);

	point_t Pt(-signed_solid_angle * Nt);
	for (unsigned int v = 0; v < 3; ++v) Pt += Nt.cross(C[v] * d[v]);
	for (unsigned int v = 0; v < 3; ++v) J[v] = e[(v + 2) % 3].cross(e[(v + 1) % 3]);

	psi = -3.0 * signed_solid_angle * signed_volume / At;
	for (unsigned int v = 0; v < 3; ++v) psi -= C[v] * J[v].dot(Nt);
	for (unsigned int v = 0; v < 3; ++v) phi[v] += Pt.dot(J[v]) / (2.0 * At);
}

Eigen::Vector4d bilinear_sheet(double u, double v)
{
	return { (1. - u) * (1. - v), u * (1. - v), u * v, (1. - u) * v };
}

template<class point_t>
point_t interpolate_on_quad(point_t* quad_vertices, double u, double v)
{
	auto const b = bilinear_sheet(u, v);
	return quad_vertices[0] * b(0) + quad_vertices[1] * b(1) + quad_vertices[2] * b(2) + quad_vertices[3] * b(3);
}

template<class point_t>
point_t quad_u_tangent(point_t* quad_vertices, double v)
{
	return (1. - v) * (quad_vertices[1] - quad_vertices[0]) + v * (quad_vertices[2] - quad_vertices[3]);
}

template<class point_t>
point_t quad_v_tangent(point_t* quad_vertices, double u)
{
	return (1. - u) * (quad_vertices[3] - quad_vertices[0]) + u * (quad_vertices[2] - quad_vertices[1]);
}


template<class point_t>
point_t interpolated_quad_normal(point_t* quad_vertices, double u, double v)
{
	auto const u_tangent = quad_u_tangent(quad_vertices, v);
	auto const v_tangent = quad_v_tangent(quad_vertices, u);
	return u_tangent.cross(v_tangent);
}

template<class point_t>
Eigen::VectorXd computePhiAndPsiForOneQuad(point_t const& eta, point_t* quad_vertices)
{
	double u_projected = 0, v_projected = 0;
	smoothProjectOnQuad(eta, quad_vertices, u_projected, v_projected);
	assert(0. <= u_projected && u_projected <= 1.);
	assert(0. <= v_projected && v_projected <= 1.);

	Eigen::VectorXd Phi;
	Phi.resize(8);
	Phi.fill(0);
	const int n = 0;
	const double m = 3.;
	const int N = 2 * n + 2;
	double u[N + 1], v[N + 1];

	for (int i = 0; i < N + 1; ++i)
	{
		if (i < N / 2)
		{
			const double factor = 1. - std::pow((static_cast<double>(N / 2) - static_cast<double>(i)) / static_cast<double>(N / 2), m);
			u[i] = factor * u_projected;
			v[i] = factor * v_projected;
		}
		else if (i == N / 2)
		{
			u[i] = u_projected;
			v[i] = v_projected;
		}
		else
		{
			const double factor = 1. - std::pow((static_cast<double>(i) - static_cast<double>(N / 2)) / static_cast<double>(N / 2), m);
			u[i] = (u_projected - 1.) * factor + 1.;
			v[i] = (v_projected - 1.) * factor + 1.;
		}
	}

	auto Proc = [&](Eigen::Vector3i u_idx, Eigen::Vector3i v_idx)
	{
		double u_tri[3], v_tri[3];
		for (int j = 0; j < 3; ++j)
		{
			u_tri[j] = u[u_idx(j)];
			v_tri[j] = v[v_idx(j)];
		}
		auto const u_avg = (u_tri[0] + u_tri[1] + u_tri[2]) / 3.;
		auto const v_avg = (v_tri[0] + v_tri[1] + v_tri[2]) / 3.;
		point_t tesselated_tri[3];
		for (int k = 0; k < 3; ++k)
		{
			tesselated_tri[k] = interpolate_on_quad(quad_vertices, u_tri[k], v_tri[k]);
		}
		auto const b_avg = bilinear_sheet(u_avg, v_avg);

		point_t Nt = (tesselated_tri[1] - tesselated_tri[0]).cross(tesselated_tri[2] - tesselated_tri[0]);
		double NtNorm = Nt.norm();
		double At = NtNorm / 2.0;
		Nt /= NtNorm;

		if (At < 0.0000000001)
		{
			return;
		}

		double psi_tri = 0.0;

		point_t e[3];    double e_norm[3];   point_t e_normalized[3];    double R[3];    point_t d[3];    double d_norm[3];     double C[3];     point_t J[3];
		for (unsigned int v = 0; v < 3u; ++v) e[v] = tesselated_tri[v] - eta;
		for (unsigned int v = 0; v < 3u; ++v) e_norm[v] = e[v].norm();
		for (unsigned int v = 0; v < 3u; ++v) e_normalized[v] = e[v] / e_norm[v];

		auto const omega_tri = get_signed_solid_angle(e_normalized[0], e_normalized[1], e_normalized[2]);
		auto const signed_solid_angle = omega_tri / (4.f * M_PI);
		auto const signed_volume = (e[0].cross(e[1])).dot(e[2]) / 6.0;

		for (unsigned int v = 0; v < 3; ++v) R[v] = e_norm[(v + 1) % 3] + e_norm[(v + 2) % 3];
		for (unsigned int v = 0; v < 3; ++v) d[v] = tesselated_tri[(v + 1) % 3] - tesselated_tri[(v + 2) % 3];
		for (unsigned int v = 0; v < 3; ++v) d_norm[v] = d[v].norm();
		for (unsigned int v = 0; v < 3; ++v) C[v] = std::log((R[v] + d_norm[v]) / (R[v] - d_norm[v])) / (4.0 * M_PI * d_norm[v]);

		point_t Pt(-signed_solid_angle * Nt);
		for (unsigned int v = 0; v < 3; ++v) Pt += Nt.cross(C[v] * d[v]);
		for (unsigned int v = 0; v < 3; ++v) J[v] = e[(v + 2) % 3].cross(e[(v + 1) % 3]);

		psi_tri = -3.0 * signed_solid_angle * signed_volume / At;
		for (unsigned int v = 0; v < 3; ++v) psi_tri -= C[v] * J[v].dot(Nt);

		auto const phi_quad = b_avg * omega_tri / (4.f * M_PI);
		auto const interpolated_normal = interpolated_quad_normal(quad_vertices, u_avg, v_avg);
		auto const psi_quad = psi_tri * b_avg / interpolated_normal.norm();
		for (unsigned int k = 0; k < 4u; ++k)
		{
			assert(phi_quad(k) == phi_quad(k));
			auto const phi_k = phi_quad(k);
			Phi(k) += phi_k;
			assert(psi_quad(k) == psi_quad(k));
			Phi(k + 4u) += psi_quad(k);
		}
	};

	for (int v_it = 0; v_it <= n; ++v_it)
	{
		for (int u_it = v_it; u_it <= (2 * n - v_it); ++u_it)
		{
			if ((u_it + v_it) % 2 == 0)
			{
				Proc({ u_it, u_it + 2, u_it + 1 }, { v_it, v_it, v_it + 1 });
				Proc({ u_it, u_it + 1, u_it + 2 }, { N - v_it, N - v_it - 1, N - v_it });
				Proc({ v_it, v_it + 1, v_it }, { u_it, u_it + 1, u_it + 2 });
				Proc({ N - v_it, N - v_it, N - v_it - 1 }, { u_it, u_it + 2, u_it + 1 });
			}
			else
			{
				Proc({ u_it, u_it + 1, u_it + 2 }, { v_it + 1, v_it, v_it + 1 });
				Proc({ u_it, u_it + 2, u_it + 1 }, { N - v_it - 1, N - v_it - 1, N - v_it });
				Proc({ v_it + 1, v_it + 1, v_it }, { u_it, u_it + 2, u_it + 1 });
				Proc({ N - v_it - 1, N - v_it, N - v_it - 1 }, { u_it, u_it + 1, u_it + 2 });
			}
		}
	}

	return Phi;
}

template<class point_t>
float_t sigma_L(point_t* old_quad_vertices, point_t* new_quad_vertices, double u, double v)
{
	auto const old_u_tangent = quad_u_tangent(old_quad_vertices, v);
	auto const old_v_tangent = quad_v_tangent(old_quad_vertices, u);
	auto const new_u_tangent = quad_u_tangent(new_quad_vertices, v);
	auto const new_v_tangent = quad_v_tangent(new_quad_vertices, u);

	return std::sqrt((new_u_tangent.squaredNorm() * old_v_tangent.squaredNorm() + old_u_tangent.squaredNorm() * new_v_tangent.squaredNorm() - 2. *
		new_u_tangent.dot(new_v_tangent) * old_u_tangent.dot(old_v_tangent)) / (2. * (old_u_tangent.cross(old_v_tangent)).squaredNorm()));
}

template<class point_t>
float_t sigma_a(point_t* old_quad_vertices, point_t* new_quad_vertices, double u, double v)
{
	auto const old_normal = interpolated_quad_normal(old_quad_vertices, u, v);
	auto const new_normal = interpolated_quad_normal(new_quad_vertices, u, v);
	return new_normal.norm() / old_normal.norm();
}

template<class point_t>
Eigen::Vector4d numerically_approx_sigma_q(point_t* old_quad_vertices, point_t* new_quad_vertices, int n = 10)
{
	auto const numSamples = n * n;
	auto const dx = 1. / static_cast<double>(n);
	Eigen::Vector4d res_sigma_a = { 0, 0, 0, 0 };
	Eigen::Vector4d res_sigma_L = { 0, 0, 0, 0 };
	for (int i = 1; i <= n; ++i)
	{
		auto const u = static_cast<double>(i) * dx - dx * .5;

		for (int j = 1; j <= n; ++j)
		{
			auto const v = static_cast<double>(j) * dx - dx * .5;

			auto const b = bilinear_sheet(u, v);
			auto const normal = interpolated_quad_normal(old_quad_vertices, u, v);
			auto const sigma_a_approx = sigma_a(old_quad_vertices, new_quad_vertices, u, v);
			auto const sigma_L_approx = sigma_L(old_quad_vertices, new_quad_vertices, u, v);
			res_sigma_a += dx * dx * b * sigma_a_approx / normal.stableNorm();
			res_sigma_L += dx * dx * b * sigma_L_approx / normal.stableNorm();
		}
	}
	Eigen::Vector4d sigma = { 0, 0, 0, 0 };
	for (int i = 0; i < 4; ++i)
	{
		sigma(i) = res_sigma_L(i) / res_sigma_a(i);
	}

	return sigma;
}

void calculateGreenCoordinatesFromQMVC(const Eigen::MatrixXd& C, const Eigen::MatrixXi& CF, const Eigen::MatrixXd& normals, const Eigen::MatrixXd& eta_m,
	Eigen::MatrixXd& phi, Eigen::MatrixXd& psi)
{
	phi.resize(C.rows(), eta_m.rows());
	psi.resize(CF.rows(), eta_m.rows());
	phi.fill(0); psi.fill(0);

	for (int i = 0; i < eta_m.rows(); ++i)
	{
		const Eigen::Vector3d eta = eta_m.row(i);

		for (unsigned int face_idx = 0; face_idx < CF.rows(); ++face_idx)
		{
			Eigen::Vector3d tri_verts[3];
			const Eigen::Vector3i tri_indices = CF.row(face_idx);
			double phi_[3];
			for (unsigned int l = 0; l < 3u; ++l)
			{
				tri_verts[l] = C.row(tri_indices[l]);
			}
			computePhiAndPsiForOneTriangle(eta, tri_verts, phi_, psi(face_idx, i));

			for (unsigned int v = 0; v < 3u; ++v)
			{
				phi(tri_indices[v], i) += phi_[v];
			}
		}
	}

	for (unsigned int i = 0; i < eta_m.rows(); ++i)
	{
		double res = 0;
		for (unsigned int j = 0; j < C.rows(); ++j)
		{
			res += phi(j, i);
		}
		assert(std::abs(1. - res) < 1e-3);
	}
}

#define EPS 0.00000000001

template <typename T> T sgn(T val) {
	return val >= EPS ? T(1) : T(-1);
}

double GCTriInt(const Eigen::Vector3d p, const Eigen::Vector3d t1, const Eigen::Vector3d t2, const Eigen::Vector3d eta = { 0., 0., 0. })
{
	auto const t2mt1 = t2 - t1;
	auto const pmt1 = p - t1;

	auto tmpVal = t2mt1.dot(p - t1) / (t2mt1.stableNorm() * pmt1.stableNorm());

	if (std::abs(tmpVal) > 1.0 - EPS)
	{
		return 0;
	}

	auto const alpha = std::acos(tmpVal);

	if (std::abs(alpha - M_PI) < EPS || std::abs(alpha) < EPS)
	{
		return 0.0;
	}

	auto const t1mp = t1 - p;
	auto const t2mp = t2 - p;

	tmpVal = t1mp.dot(t2mp) / (t1mp.stableNorm() * t2mp.stableNorm());

	if (std::abs(tmpVal) > 1.0 - EPS)
	{
		return 0;
	}

	auto const beta = std::acos(tmpVal);

	auto const lambda = pmt1.squaredNorm() * std::sin(alpha) * std::sin(alpha);
	auto const pmeta = p - eta;
	auto const c = pmeta.squaredNorm();

	auto const sqrt_c = std::sqrt(c);
	auto const sqrt_lambda = std::sqrt(lambda);

	auto getI = [&](double theta) {
		auto const S = std::sin(theta);
		auto const C = std::cos(theta);
		auto const fac_sign = -.5 * static_cast<double>(sgn(S));
		auto const add_1 = 2. * sqrt_c * std::atan(C * sqrt_c / std::sqrt(lambda + S * S * c));
		auto const denom1 = (1. - C) * (1. - C);
		auto const denom2 = c * (1. + C) + lambda + std::sqrt(lambda * lambda + lambda * c * S * S);
		auto const add_2 = sqrt_lambda * std::log((2. * sqrt_lambda * S * S) / denom1 * (1. - 2. * c * C / denom2));
		return fac_sign * (add_1 + add_2);
	};

	auto const I_1 = getI(M_PI - alpha);
	auto const I_2 = getI(M_PI - alpha - beta);

	return (-1. / (4 * M_PI)) * std::abs(I_1 - I_2 - sqrt_c * beta);
}

void calculateGreenCoordinates(const Eigen::MatrixXd& C, const Eigen::MatrixXi& CF, const Eigen::MatrixXd& normals, const Eigen::MatrixXd& eta_m,
	Eigen::MatrixXd& phi, Eigen::MatrixXd& psi)
{
	phi.resize(C.rows(), eta_m.rows());
	psi.resize(CF.rows(), eta_m.rows());
	phi.fill(0); psi.fill(0);

	for (int i = 0; i < eta_m.rows(); ++i)
	{
		const Eigen::Vector3d eta = eta_m.row(i);

		for (int face_idx = 0; face_idx < CF.rows(); ++face_idx)
		{
			const Eigen::Vector3i index_vec = CF.row(face_idx);
			const Eigen::Vector3d t_0 = C.row(index_vec[0]);
			const Eigen::Vector3d t_1 = C.row(index_vec[1]);
			const Eigen::Vector3d t_2 = C.row(index_vec[2]);

			Eigen::Vector3d v[3];

			v[0] = t_0 - eta;
			v[1] = t_1 - eta;
			v[2] = t_2 - eta;

			const Eigen::Vector3d normal = normals.row(face_idx);
			assert(std::abs(1.0 - normal.stableNorm()) < 1e-4);

			auto const p = v[0].dot(normal) * normal;

			double s[3];
			double I[3];
			double II[3];
			Eigen::Vector3d N[3];

			for (unsigned int l = 0; l < 3u; ++l)
			{
				s[l] = sgn(((v[l] - p).cross(v[(l + 1u) % 3u] - p)).dot(normal));
				I[l] = GCTriInt(p, v[l], v[(l + 1) % 3]);
				II[l] = GCTriInt({ 0., 0., 0. }, v[(l + 1) % 3], v[l]);
				if (II[l] != II[l])
				{
					std::cout << "NAN";
				}

				auto const q = v[(l + 1) % 3].cross(v[l]);
				N[l] = q.stableNormalized();
			}

			auto const I_final = -1. * std::abs(s[0] * I[0] + s[1] * I[1] + s[2] * I[2]);
			auto const psi_normal = -I_final;

			if (psi_normal != psi_normal)
			{
				std::cout << "NAN\n";
			}

			psi(face_idx, i) = psi_normal;

			Eigen::Vector3d const w = normal * I_final + N[0] * II[0] + N[1] * II[1] + N[2] * II[2];

			auto const w_norm = w.stableNorm();

			if (w.stableNorm() > 1.e-3)
			{
				for (unsigned int l = 0; l < 3u; ++l)
				{
					auto const phi_vert = (N[(l + 1u) % 3u].dot(w)) / (N[(l + 1u) % 3u].dot(v[l]));
					auto const phi_vert_abs = std::abs(phi_vert);
					phi(index_vec(l), i) += phi_vert;
				}
			}

		}
	}

	for (unsigned int i = 0; i < eta_m.rows(); ++i)
	{
		double res = 0;
		for (unsigned int j = 0; j < C.rows(); ++j)
		{
			res += phi(j, i);
		}
		assert(std::abs(1. - res) < 1e-3);
	}
}

void calculateGreenCoordinatesTriQuad(const Eigen::MatrixXd& C, const Eigen::MatrixXi& CF, Eigen::MatrixXd const& eta_m,
	Eigen::MatrixXd& phi, std::vector<double>& psi_tri, std::vector<Eigen::Vector4d>& psi_quad)
{
	phi.resize(C.rows(), eta_m.rows());
	phi.fill(0);

	bool contains_quads = CF.cols() == 4;

	for (int eta_idx = 0; eta_idx < eta_m.rows(); ++eta_idx)
	{
		const Eigen::Vector3d eta = eta_m.row(eta_idx);

		for (int face_idx = 0; face_idx < CF.rows(); ++face_idx)
		{
			auto const face = CF.row(face_idx);
			Eigen::Vector3d face_points[4];

			for (int k = 0; k < face.size(); ++k)
			{
				face_points[k] = C.row(face(k));
			}

			if (face.size() == 3 || (contains_quads && face(3) == -1))
			{
				double phi_[3];
				double psi = 0;
				computePhiAndPsiForOneTriangle(eta, face_points, phi_, psi);

				for (int k = 0; k < 3; ++k)
				{
					phi(face(k), eta_idx) += phi_[k];
				}
				psi_tri.push_back(psi);
			}
			else if (face.size() == 4)
			{
				auto const Phi = computePhiAndPsiForOneQuad(eta, face_points);
				assert(Phi.size() == 8);

				for (int k = 0; k < 4; ++k)
				{
					phi(face[k], eta_idx) += Phi(k);
				}
				psi_quad.push_back({ Phi(4), Phi(5), Phi(6), Phi(7) });
			}
			else
			{
				std::cerr << "QGC supports only triangles and quads! Unsupported face type found!\n";
			}
		}
		/*double res = 0;
		for (unsigned int j = 0; j < C.rows(); ++j)
		{
			res += phi(j, eta_idx);
		}
		assert(std::abs(1. - res) < 1e-3);*/
	}
}

void calcNormals(const Eigen::MatrixXd& C, const Eigen::MatrixXi& CF, Eigen::MatrixXd& normals)
{
	normals.resize(CF.rows(), 3);

	for (int i = 0; i < CF.rows(); ++i)
	{
		Eigen::Vector3i index_vector = CF.row(i);

		const Eigen::Vector3d t_0 = C.row(index_vector[0]);
		const Eigen::Vector3d t_1 = C.row(index_vector[1]);
		const Eigen::Vector3d t_2 = C.row(index_vector[2]);

		auto const normal = ((t_1 - t_0).cross(t_2 - t_0)).normalized();

		normals.row(i) = normal;
	}
}

double calc_scaling_factor_tri(Eigen::Vector3d old_tri[3], Eigen::Vector3d new_tri[3])
{
	auto const old_u = old_tri[0] - old_tri[1];
	auto const old_v = old_tri[0] - old_tri[2];

	auto const area = .5 * (old_u.cross(old_v)).stableNorm();

	auto const new_u = new_tri[0] - new_tri[1];
	auto const new_v = new_tri[0] - new_tri[2];

	return std::sqrt(new_u.squaredNorm() * old_v.squaredNorm() - 2. * new_u.dot(new_v) * old_u.dot(old_v) + new_v.squaredNorm() * old_u.squaredNorm()) /
		(2.8284271247461903 * area);
}

void calcScalingFactors(const Eigen::MatrixXd& C, const Eigen::MatrixXd& C_deformed, const Eigen::MatrixXi& CF, Eigen::MatrixXd& normals)
{
	for (int i = 0; i < CF.rows(); ++i)
	{
		const Eigen::Vector3i index_vec = CF.row(i);

		Eigen::Vector3d old_tri[3], new_tri[3];
		for (int k = 0; k < 3; ++k)
		{
			old_tri[k] = C.row(index_vec(k));
			new_tri[k] = C_deformed.row(index_vec(k));
		}

		auto const scaling_factor = calc_scaling_factor_tri(old_tri, new_tri);
		normals.row(i) *= scaling_factor;
	}
}

void calcNewPositionsTriQuad(const Eigen::MatrixXd& C, const Eigen::MatrixXd& C_deformed, const Eigen::MatrixXi & CF,
	Eigen::MatrixXd& phi, std::vector<double> const& psi_tri, std::vector<Eigen::Vector4d> const& psi_quad, Eigen::MatrixXd& eta_deformed)
{
	// Adding weighted control vertices first
	eta_deformed = phi.transpose() * C_deformed;

	// Precalculate scaling factors of quads
	std::vector<Eigen::Vector4d> sigma_quads;
	for (unsigned int face_idx = 0; face_idx < CF.rows(); ++face_idx)
	{
		auto const face = CF.row(face_idx);

		if (face.size() != 4 || face(3) != -1)
		{
			continue;
		}
		Eigen::Vector3d new_face_points[4], old_face_points[4];
		for (int k = 0; k < face.size(); ++k)
		{
			old_face_points[k] = C.row(face(k));
			new_face_points[k] = C_deformed.row(face(k));
		}

		auto const sigma_q = numerically_approx_sigma_q(old_face_points, new_face_points);
		sigma_quads.push_back(sigma_q);
	}

	// Now adding psi for every face
	unsigned int tri_count = 0, quad_count = 0;
	for (unsigned int eta_idx = 0; eta_idx < eta_deformed.rows(); ++eta_idx)
	{
		Eigen::Vector3d eta = eta_deformed.row(eta_idx);
		unsigned int loc_quad_count = 0;

		for (unsigned int face_idx = 0; face_idx < CF.rows(); ++face_idx)
		{
			auto const face = CF.row(face_idx);
			Eigen::Vector3d new_face_points[4], old_face_points[4];

			for (int k = 0; k < face.size(); ++k)
			{
				old_face_points[k] = C.row(face(k));
				new_face_points[k] = C_deformed.row(face(k));
			}

			if (face.size() == 3 || (face.size() == 4 && face(3) == -1))
			{
				auto const sigma_t = calc_scaling_factor_tri(old_face_points, new_face_points);
				auto const psi_t = psi_tri[tri_count++];
				auto const new_normal = ((new_face_points[1] - new_face_points[0]).cross(new_face_points[2] - new_face_points[0])).normalized();
				eta += sigma_t * psi_t * new_normal;
			}
			else if (face.size() == 4)
			{
				auto const sigma_q = sigma_quads[loc_quad_count++];
				auto const psi_q = psi_quad[quad_count++];
				for (int k = 0; k < 4; ++k)
				{
					auto const new_N_k = (new_face_points[(k + 1) % 4] - new_face_points[k]).cross(new_face_points[(k + 3) % 4] - new_face_points[k]);
					eta += psi_q(k) * sigma_q(k) * new_N_k;
				}
			}
			else
			{
				std::cerr << "QGC supports only triangles and quads! Unsupported face type found!\n";
			}

		}

		eta_deformed.row(eta_idx) = eta;
	}
}