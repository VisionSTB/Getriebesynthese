#ifndef __GRUNDPUNKT_H
#define __GRUNDPUNKT_H

/*
#include <fstream>
#include <iostream>     // std::cout
#include <algorithm>
#include <regex>
#include <map>
#include <utility>      // std::pair

#include <string>       // std::string
#include <streambuf>
#include <sstream>      // std::ostringstream
*/
#include <vector>
#include <array>

#include <SDL.h>

//using namespace std::string_literals;

struct SPointD
    {
    SPointD() = default;
    template<typename T>
	SPointD(T x, T y) : x(x), y(y) {}
    SPointD(SDL_Point const & P) : x(P.x), y(P.y) {}

    double x{0}, y{0};
    explicit operator SDL_Point()
	{
	return { (int)x, (int)y };
	}
    SPointD operator - ( SPointD const & P ) const { return { x-P.x, y-P.y }; }
    };

struct SLineD
    {
    bool   bHasP1{false};
    double x1{0}, y1{0};
    double x2{0}, y2{0};

    SLineD operator - ( SPointD const & P ) const { return { bHasP1, x1-P.x, y1-P.y, x2-P.x, y2-P.y }; }
    };



struct SUmkreisD
    {
    double  R;
    SPointD M;
    };

struct SCollision
    {
    int E; // 1,2,3 = E1, E2, E3
    int P; // 1,2,3 = P1, P2, Pm
    };


//using VLines       = std::vector<SLineD>;
using VEbenenLagen = std::vector<SLineD>;
using VPolDreieck  = std::vector<SPointD>;
using VGelenke     = std::vector<SPointD>;
using A3Gelenke    = std::array<SPointD, 3>;

SLineD FixedLenLine(SLineD & roL, double const & crnLenEbene, bool const & crbFirst = true)
    {
    double const dx   = roL.x1 - roL.x2;
    double const dy   = roL.y1 - roL.y2;
    double const nLen = sqrt(dx*dx + dy*dy);
    double const q    = (double)crnLenEbene / ((nLen!=0)?nLen:1);
    if (crbFirst)
        {
        roL.x2 = roL.x1 - dx*q;
        roL.y2 = roL.y1 - dy*q;
        }
    else
        {
        roL.x1 = roL.x2 + dx*q;
        roL.y1 = roL.y2 + dy*q;
        }
    return roL;
    } // void FixedLenLine(...

SPointD Intersection(SLineD const & E1, SLineD const & E2)
    {
    double dx1,dx2,m1,n1,m2,n2;

    dx1 = E1.x2 - E1.x1;
    dx2 = E2.x2 - E2.x1;

    m1 = (E1.y2 - E1.y1) / dx1; // Steigungen ermitteln
    m2 = (E2.y2 - E2.y1) / dx2;

    // if (ROUND(m1,MAX_ACCURACY)==ROUND(m2,MAX_ACCURACY)) return false; // Geraden sind parallel

    n1 = E1.y1 - (m1*E1.x1); // Abstände von X-Achse ermitteln
    n2 = E2.y1 - (m2*E2.x1);

    double x = (n2-n1)/(m1-m2); // Schnittpunktkoordinate berechnen
    double y = m1*x+n1;

    return { x, y };
  } // Intersection


SLineD Perpendicle(SLineD const & croLine)
    {
    auto const dx = (croLine.x2 - croLine.x1)/2.0;
    auto const dy = (croLine.y2 - croLine.y1)/2.0;

    SLineD I;
    I.x1 = croLine.x2 - dy - dx;
    I.y1 = croLine.y2 + dx - dy;
    I.x2 = croLine.x2 + dy - dx;
    I.y2 = croLine.y2 - dx - dy;

    return std::move(I);
    }

SPointD PointMirror(SDL_Renderer * pSdlRenderer, SPointD const & croPoint, SLineD const & croMirror)
    {
    auto const dx = (croMirror.x2 - croMirror.x1)/2.0;
    auto const dy = (croMirror.y2 - croMirror.y1)/2.0;

    SLineD I;
    I.x1 = croPoint.x - dy - 0*dx;
    I.y1 = croPoint.y + dx - 0*dy;
    I.x2 = croPoint.x + dy - 0*dx;
    I.y2 = croPoint.y - dx - 0*dy;
//  SDL_RenderDrawLine( pSdlRenderer,  I.x1,  I.y1,  I.x2,  I.y2 );

    auto const S  = Intersection(croMirror, I);
    auto const mx = (croMirror.x2 + croMirror.x1)/2.0;
    auto const my = (croMirror.y2 + croMirror.y1)/2.0;

    return { S.x + (S.x - croPoint.x), S.y + (S.y - croPoint.y) };
    }

SPointD CalcPolpunkt(SLineD const & E1, SLineD const & E2)
    {
    SLineD const L1{ true, E1.x1, E1.y1, E2.x1, E2.y1 };
    auto La = Perpendicle(L1);

    SLineD const L2{ true, E1.x2, E1.y2, E2.x2, E2.y2 };
    auto Lb = Perpendicle(L2);

    return Intersection( La, Lb );
    }


long   const g_cnSteps = 180;
using        K180      = SDL_Point[g_cnSteps];
K180         g_cnCircle;
double const g_nf{1000.0};

void MakeCircle()
    {
    long double const pi = 4*atan(1);
    for ( int nStep = 0; nStep < g_cnSteps; ++nStep )
        {
        g_cnCircle[nStep]={ (int)(cos((double)nStep/(double)(g_cnSteps-1)*2.0*pi)*g_nf), (int)(sin((double)nStep/(double)(g_cnSteps-1)*2.0*pi)*g_nf) };
        }
    } // void MakeCircle()

void RenderCircle(SDL_Renderer * pSdlRenderer, int const & crnRadius, int const & X, int const & Y )
    {
    K180 aKreis;
    memcpy(aKreis, g_cnCircle, sizeof(g_cnCircle));
    for ( int step = 0; step < g_cnSteps; ++step )
        {
        aKreis[step]={(int)((double)(aKreis[step].x)*(double)crnRadius/g_nf+X), (int)((double)(aKreis[step].y)*(double)crnRadius/g_nf+Y)};
        }
    SDL_RenderDrawLines(pSdlRenderer, aKreis, g_cnSteps);
    } // void RenderCircle(...

SUmkreisD RenderUmkreis( SDL_Renderer * pSdlRenderer, SPointD const & P1, SPointD const & P2, SPointD const & P3 )
    {
    SLineD const L1{ true, (double)P1.x, (double)P1.y, (double)P2.x, (double)P2.y };
    auto La = Perpendicle(L1);
    SLineD const L2{ true, (double)P1.x, (double)P1.y, (double)P3.x, (double)P3.y };
    auto Lb = Perpendicle(L2);

    auto M = Intersection( La, Lb );
    auto R = sqrt( (double)(M.x - L1.x1)*(M.x - L1.x1) + (M.y - L1.y1)*(M.y - L1.y1) );

    RenderCircle(pSdlRenderer, R, M.x, M.y );

    return { R, M };
    }






class CGrundpunkt
    {
    protected:

	double    m_dX{};
	double    m_dY{};
	SUmkreisD m_tUK{};
	A3Gelenke m_a3Gelenke{}; // Gelenkpunkte
	bool      m_bFixed{true};


    public:

	CGrundpunkt(SPointD const & P123) : m_dX(P123.x), m_dY(P123.y) {}
	CGrundpunkt(CGrundpunkt const & src) = default;
/*
	CGrundpunkt(CGrundpunkt const & src)
	    {
	    m_dX        = src.m_dX;
	    m_dY        = src.m_dY;
	    m_tUK.R       = src.m_tUK.R;
	    m_tUK.M       = src.m_tUK.M;
	    m_a3Gelenke = src.m_a3Gelenke;
	    m_bFixed    = src.m_bFixed;
	    }
*/
	void FixIt( bool bFixit = true ) { m_bFixed = bFixit; }
	bool Isfix() { return m_bFixed; }

	SPointD const & G0() const { return m_tUK.M; }
	SPointD const & GPoint( int i ) const { return m_a3Gelenke[i]; }

	void UpdateAndShow(SDL_Renderer * pSdlRenderer, SPointD const & P123, VPolDreieck const & Poldreieck)
	    {
	    Update(P123);
	    Show(pSdlRenderer, Poldreieck);
	    }

	void Update(SPointD const & P123)
	    {
	    m_dX = P123.x;
	    m_dY = P123.y;
	    }

	void Show(SDL_Renderer * pSdlRenderer, VPolDreieck const & Poldreieck)
	    {
	    RenderCircle( pSdlRenderer, 8, m_dX, m_dY );
	    RenderCircle( pSdlRenderer, 7, m_dX, m_dY );
	    RenderCircle( pSdlRenderer, 6, m_dX, m_dY );
	    RenderCircle( pSdlRenderer, 5, m_dX, m_dY );

	    A3Gelenke G;
	    for ( int n=0, i=0, j=0; n < 3; ++n )
		{
		// Lines: P12-P13, P12-P23, P13-P23
		if (n == 0) { i=0; j=1; }
		if (n == 1) { i=0; j=2; }
		if (n == 2) { i=1; j=2; }
		SLineD PL{ true, (double)Poldreieck[i].x, (double)Poldreieck[i].y,
				 (double)Poldreieck[j].x, (double)Poldreieck[j].y};
		G[n] = PointMirror( pSdlRenderer, { m_dX, m_dY }, PL );
		RenderCircle( pSdlRenderer, 20, G[n].x, G[n].y );
		} // for (...
	    m_a3Gelenke = G;

	    m_tUK = RenderUmkreis(pSdlRenderer, G[0], G[1], G[2] );
	    RenderCircle( pSdlRenderer, 36, m_tUK.M.x, m_tUK.M.y );

	    for (auto const & a:G)
		{
		SDL_RenderDrawLine(pSdlRenderer, m_tUK.M.x, m_tUK.M.y, a.x, a.y);
		}
	    }
    }; // class CGrundpunkt

using VGrundpunkte  = std::vector<CGrundpunkt>;
using A2Grundpunkte = std::array<CGrundpunkt, 2>;

// __GRUNDPUNKT_H
#endif
