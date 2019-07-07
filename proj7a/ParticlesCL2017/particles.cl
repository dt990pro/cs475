typedef float4 point;
typedef float4 vector;
typedef float4 color;
typedef float4 sphere;


vector
Bounce( vector in, vector n )
{
	vector out = in - n*(vector)( 2.*dot(in.xyz, n.xyz) );
	out.w = 0.;
	return out;
}

vector
BounceSphere( point p, vector v, sphere s )
{
	vector n;
	n.xyz = fast_normalize( p.xyz - s.xyz );
	n.w = 0.;
	return Bounce( v, n );
}

bool
IsInsideSphere( point p, sphere s )
{
	float r = fast_length( p.xyz - s.xyz );
	return  ( r < s.w );
}

kernel
void
Particle( global point *dPobj, global vector *dVel, global color *dCobj )
{
	float4 G;
	const float  DT      = 0.1;
	const sphere Sphere1 = (sphere)( -350., 0., 0.,  300. );
	const sphere Sphere2 = (sphere)( 350., 0., 0.,  300. );

	int gid = get_global_id( 0 );

	// ?????
	point p = dPobj[gid];
	vector v = dVel[gid];

	// update 1 DT
	if (p.y > 0.){
		G       = (float4) ( 0., -17.6, 0., 0. );
	}
	else{
		G       = (float4) ( 0., 17.6, 0., 0. );
	}
	point pp = p + v*DT + G*(point)(0.5*DT*DT);
	vector vp = v + G*DT;
	pp.w = 1.;
	vp.w = 0.;

	if (IsInsideSphere(pp, Sphere1)){
		// update 1 DT
		vp = BounceSphere(p, v, Sphere1);
		pp = p + vp*DT + G*(point)(0.5*DT*DT);

		dCobj[gid] = (vector)(.9, .6, .3, 1.);
	}

	if (IsInsideSphere(pp, Sphere2)){
		// update 1 DT
		vp = BounceSphere(p, v, Sphere2);
		pp = p + vp*DT + G*(point)(0.5*DT*DT);

		dCobj[gid] = (vector)(.3, .6, .9, 1.);
	}

	dPobj[gid] = pp;
	dVel[gid] = vp;
}
