#include <string>
#include <iostream>
#include <vector>
#include <array>
#include <functional> // std::optional

#include <SDL.h>

#include "grundpunkt.h"
#include "te.h"


struct SLineD
    {
    bool   bHasP1{false};
    double x1{0}, y1{0};
    double x2{0}, y2{0};
    };

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
    };

struct SUmkreisD
    {
    double  R;
    SPointD M;
    };

//using VLines       = std::vector<SLineD>;
using VEbenenLagen = std::vector<SLineD>;
using VPolDreieck  = std::vector<SPointD>;
using VGelenke     = std::vector<SPointD>;
using A3Gelenke    = std::array<SPointD, 3>;


//VLines       g_vLines;
VEbenenLagen g_vEbenen;     // E1,  E2,  E3  = 3 homologe Lagen einer Ebene
VPolDreieck  g_vPoldreieck; // P12, P13, P23 = 3 Polpunkte zu den o.g.Ebenenlagen
VGelenke     g_vGelenke;

struct SCollision
    {
    int E; // 1,2,3 = E1, E2, E3
    int P; // 1,2,3 = P1, P2, Pm
    } g_tCollision;

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
	SUmkreisD m_tUK;
	A3Gelenke m_a3Gelenke; // Gelenkpunkte
	bool      m_bFixed{true};


    public:

	CGrundpunkt(SPointD const & P123) : m_dX(P123.x), m_dY(P123.y) {}

	void FixIt( bool bFixit = true ) { m_bFixed = bFixit; }
	bool Isfix() { return m_bFixed; }

	SPointD GPoint( int i ) { return m_a3Gelenke[i]; }

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

	    SUmkreisD m_tUK = RenderUmkreis(pSdlRenderer, G[0], G[1], G[2] );
	    RenderCircle( pSdlRenderer, 36, m_tUK.M.x, m_tUK.M.y );

	    for (auto const & a:G)
		{
		SDL_RenderDrawLine(pSdlRenderer, m_tUK.M.x, m_tUK.M.y, a.x, a.y);
		}
	    }
/*
	if ( g_vGelenke.size() < 3 )
	    {
	    g_vGelenke.emplace_back( SDL_Point{(int)G12.x, (int)G12.y} );
	    g_vGelenke.emplace_back( SDL_Point{(int)G13.x, (int)G13.y} );
	    g_vGelenke.emplace_back( SDL_Point{(int)G23.x, (int)G23.y} );
	    }
	else
	    {
	    g_vGelenke[0] = {(int)G12.x, (int)G12.y};
	    g_vGelenke[1] = {(int)G13.x, (int)G13.y};
	    g_vGelenke[2] = {(int)G23.x, (int)G23.y};
	    }
*/
    }; // class CGrundpunkt

CGrundpunkt g_oGrundpunkt ( {-1, -1} );

using VGrundpunkte  = std::vector<CGrundpunkt>;
using A2Grundpunkte = std::array<CGrundpunkt, 2>;

VGrundpunkte  g_vGrundpunkte;  // A123, B123


/******************************************************************************************

*/

int main(int argc, char* argv[])
    {
    MakeCircle();

    if (SDL_Init(SDL_INIT_VIDEO) == 0)
        {
        SDL_Window *   pSdlWindow   = nullptr;
        SDL_Renderer * pSdlRenderer = nullptr;

        if (SDL_CreateWindowAndRenderer(800, 600, 0, &pSdlWindow, &pSdlRenderer) == 0)
            {
            SDL_SetWindowFullscreen(pSdlWindow, SDL_WINDOW_FULLSCREEN);
            SDL_bool done = SDL_FALSE;

            SLineD  oL;
            bool    bLeftMouseDown{true};
            double  nLenEbene{0};
            while (!done)
                {
                SDL_PumpEvents();
                int x, y;
                if ( SDL_GetMouseState(&x, &y) & SDL_BUTTON(SDL_BUTTON_LEFT) )
                    {
                    bLeftMouseDown = false;
                    if ( oL.bHasP1 )
                        {
                        oL.x2 = x;
                        oL.y2 = y;
                        }
                    else
                        {
                        oL.x1 = x;
                        oL.y1 = y;
                        oL.bHasP1 = true;
                        }
                    }
                else
                    if ( !bLeftMouseDown )
                        {

                        if (!g_oGrundpunkt.Isfix())
			    {
			    g_oGrundpunkt.FixIt();
			    if ( g_vGrundpunkte.size() < 2 )
				{
				g_vGrundpunkte.emplace_back(g_oGrundpunkt);
				}
			    }

                	if ( (g_vEbenen.size() < 3) )
                            {
                            if ( g_vEbenen.size() > 0 ) FixedLenLine(oL, nLenEbene);
                            g_vEbenen.emplace_back(oL);
                            switch ( g_vEbenen.size() )
                                {
                                case 2: g_vPoldreieck.emplace_back( CalcPolpunkt(g_vEbenen[1-1], g_vEbenen[2-1]) );
                                    break;
                                case 3: g_vPoldreieck.emplace_back( CalcPolpunkt(g_vEbenen[1-1], g_vEbenen[3-1]) );
                                        g_vPoldreieck.emplace_back( CalcPolpunkt(g_vEbenen[2-1], g_vEbenen[3-1]) );
                                     break;
                                }
                            }
                        else
                            {
                        //-    g_vLines.emplace_back(oL);
                            }
                        oL.bHasP1 = false;
                        bLeftMouseDown = true;
                        }

                SDL_SetRenderDrawColor(pSdlRenderer, 255, 255, 255, SDL_ALPHA_OPAQUE);
                SDL_RenderClear(pSdlRenderer);

                SDL_SetRenderDrawColor(pSdlRenderer, 255, 0, 0, SDL_ALPHA_OPAQUE);
                if ( !bLeftMouseDown )
                    {
                    if ( (g_vEbenen.size() > 2) || (g_vEbenen.size() < 1) )
                        {
                        if (g_vEbenen.size() < 1) SDL_RenderDrawLine(pSdlRenderer, oL.x1, oL.y1, oL.x2, oL.y2);

                        // post-move
                        if ( (g_tCollision.E > 0) && (g_tCollision.P > 0) )
                            {
                            switch ( g_tCollision.P )
                                {
                                case 1: g_vEbenen[g_tCollision.E-1].x1 = x; g_vEbenen[g_tCollision.E-1].y1 = y; break;
                                case 2: g_vEbenen[g_tCollision.E-1].x2 = x; g_vEbenen[g_tCollision.E-1].y2 = y; break;
                                }
                            SDL_RenderDrawLine(pSdlRenderer, g_vEbenen[g_tCollision.E-1].x1, g_vEbenen[g_tCollision.E-1].y1, g_vEbenen[g_tCollision.E-1].x2, g_vEbenen[g_tCollision.E-1].y2);
                            }
                        else
                            {
                            if ( g_vPoldreieck.size() > 2 )
                        	{
                        	g_oGrundpunkt.FixIt(false);
                        	g_oGrundpunkt.UpdateAndShow(pSdlRenderer, {(double)x, (double)y}, g_vPoldreieck);

            		    if ( g_vGrundpunkte.size() > 0 )
            			{
            			for ( int i = 0; i < 3; ++i )
            			    {
            			    SDL_RenderDrawLine( pSdlRenderer,
            						g_vGrundpunkte[0].GPoint(i).x, g_vGrundpunkte[0].GPoint(i).y,
            						g_oGrundpunkt.GPoint(i).x,    g_oGrundpunkt.GPoint(i).y );
            			    }
            			}

                        	}
                            }
                        if (g_vEbenen.size() > 2)
                            {
                            FixedLenLine(g_vEbenen[0], nLenEbene, g_tCollision.P == 1);
                            FixedLenLine(g_vEbenen[1], nLenEbene, g_tCollision.P == 1);
                            FixedLenLine(g_vEbenen[2], nLenEbene, g_tCollision.P == 1);
                            g_vPoldreieck[0] = ( CalcPolpunkt(g_vEbenen[1-1], g_vEbenen[2-1]) );
                            g_vPoldreieck[1] = ( CalcPolpunkt(g_vEbenen[1-1], g_vEbenen[3-1]) );
                            g_vPoldreieck[2] = ( CalcPolpunkt(g_vEbenen[2-1], g_vEbenen[3-1]) );
                            }
                        }
                    else
                        {
                        FixedLenLine(oL, nLenEbene);
                        SDL_RenderDrawLine(pSdlRenderer, oL.x1, oL.y1, oL.x2, oL.y2);
                        }
                    }

                /// doit
                {
                SPointD tP12, tP23;
                SDL_Rect  tRect;
                switch ( g_vEbenen.size() )
                    {
                    case 1: tP12 = CalcPolpunkt(g_vEbenen[1-1], oL);
                            tRect={(int)tP12.x-8, (int)tP12.y-8, 16, 16};
                            SDL_RenderDrawRect(pSdlRenderer, &tRect);

                            break;

                    case 2: tP12 = CalcPolpunkt(g_vEbenen[1-1], oL);
                            tRect={(int)tP12.x-8, (int)tP12.y-8, 16, 16};
                            SDL_RenderDrawRect(pSdlRenderer, &tRect);
                            tP23 = CalcPolpunkt(g_vEbenen[2-1], oL);
                            tRect={(int)tP23.x-8, (int)tP23.y-8, 16, 16};
                            SDL_RenderDrawRect(pSdlRenderer, &tRect);

                            SDL_PumpEvents();
                            if ( SDL_GetMouseState(&x, &y) & SDL_BUTTON(SDL_BUTTON_LEFT) )
                                {
                                SDL_RenderDrawLine(pSdlRenderer, (int)g_vPoldreieck[0].x, (int)g_vPoldreieck[0].y, (int)tP12.x, (int)tP12.y);
                                SDL_RenderDrawLine(pSdlRenderer, (int)tP12.x, (int)tP12.y,                               (int)tP23.x, (int)tP23.y);
                                SDL_RenderDrawLine(pSdlRenderer, (int)tP23.x, (int)tP23.y, (int)g_vPoldreieck[0].x, (int)g_vPoldreieck[0].y);

                                RenderUmkreis( pSdlRenderer, g_vPoldreieck[0], tP12, tP23 );
                                }
                            break;
                    }
                SDL_SetRenderDrawColor(pSdlRenderer, 127, 127, 127, SDL_ALPHA_OPAQUE);

        	if ( g_vGrundpunkte.size() < 2 ) // ----------------------------------VIS
                if ( g_vPoldreieck.size() > 2 )
                    {
                    SDL_RenderDrawLine(pSdlRenderer, g_vPoldreieck[0].x, g_vPoldreieck[0].y, g_vPoldreieck[1].x, g_vPoldreieck[1].y);
                    SDL_RenderDrawLine(pSdlRenderer, g_vPoldreieck[1].x, g_vPoldreieck[1].y, g_vPoldreieck[2].x, g_vPoldreieck[2].y);
                    SDL_RenderDrawLine(pSdlRenderer, g_vPoldreieck[2].x, g_vPoldreieck[2].y, g_vPoldreieck[0].x, g_vPoldreieck[0].y);

                    RenderUmkreis( pSdlRenderer, g_vPoldreieck[0], g_vPoldreieck[1], g_vPoldreieck[2] );
                    }
                }

        	if ( g_vGrundpunkte.size() < 2 ) // ----------------------------------VIS
                for (auto const & a:g_vPoldreieck)
                    {
                    SDL_Rect tRect{(int)a.x-8, (int)a.y-8, 16, 16};
                    SDL_RenderDrawRect(pSdlRenderer, &tRect);
                    }

                SDL_SetRenderDrawColor(pSdlRenderer, 0, 0, 255, SDL_ALPHA_OPAQUE);
        	if ( g_vGrundpunkte.size() < 2 ) // ----------------------------------VIS
                if ( (g_vEbenen.size() > 0) && (nLenEbene == 0) )
                    {
                    int dx = g_vEbenen[0].x1 - g_vEbenen[0].x2;
                    int dy = g_vEbenen[0].y1 - g_vEbenen[0].y2;
                    nLenEbene = sqrt( dx*dx + dy*dy );
                    }

                for (auto const & a:g_vEbenen)
                    {
                    SDL_RenderDrawLine(pSdlRenderer, a.x1, a.y1, a.x2, a.y2);

                    SDL_Rect tRect{(int)a.x1-6, (int)a.y1-6, 12, 12};
                    SDL_RenderDrawRect(pSdlRenderer, &tRect);
                             tRect = {(int)a.x2-6, (int)a.y2-6, 12, 12};
                    SDL_RenderDrawRect(pSdlRenderer, &tRect);
                    }

                if ( g_vPoldreieck.size() > 2 )
                    {
		    for (auto & a:g_vGrundpunkte)
			{
			a.Show(pSdlRenderer, g_vPoldreieck);
			}
		    if ( g_vGrundpunkte.size() > 1 )
			{
			for ( int i = 0; i < 3; ++i )
			    {
			    SDL_RenderDrawLine( pSdlRenderer,
						g_vGrundpunkte[0].GPoint(i).x, g_vGrundpunkte[0].GPoint(i).y,
						g_vGrundpunkte[1].GPoint(i).x, g_vGrundpunkte[1].GPoint(i).y );
			    }
			}
                    }
/*
                SDL_SetRenderDrawColor(pSdlRenderer, 0, 0, 0, SDL_ALPHA_OPAQUE);
                for (auto const & a:g_vLines)
                    {
                    SDL_RenderDrawLine(pSdlRenderer, a.x1, a.y1, a.x2, a.y2);
                    }
*/
                SDL_SetRenderDrawColor(pSdlRenderer, 127, 0, 0, SDL_ALPHA_OPAQUE);
		if ( g_vGelenke.size() >= 3 )
		    {
		    for (auto const & a:g_vGelenke)
			RenderCircle( pSdlRenderer, 20, a.x, a.y );
		    }


                {
                double dM = 20;
                double xd = 1e6;
                double yd = 1e6;
                int     c = 0;
                int     i = 0;  // E1,2,3
                int     j = 0;  // P1,2,m
                for (auto const & a:g_vEbenen)
                    {
                    c++;
                    auto dN = sqrt( (a.x1 - x)*(a.x1 - x) + (a.y1 - y)*(a.y1 - y) );
                    if ( dN < dM )
                        {
                        dM = dN; xd = a.x1; yd = a.y1; i = c; j = 1;
                        }
                    dN = sqrt( (a.x2 - x)*(a.x2 - x) + (a.y2 - y)*(a.y2 - y) );
                    if ( dN < dM )
                        {
                        dM = dN; xd = a.x2; yd = a.y2; i = c; j = 2;
                        }
/*
                    dN = sqrt( (x-(a.x2 + a.x1)/2)*(x-(a.x2 + a.x1)/2) + (y-(a.y2 + a.y1)/2)*(y-(a.y2 + a.y1)/2) );
                    if ( dN < dM )
                        {
                        dM = dN; xd = (a.x2 + a.x1)/2; yd = (a.y2 + a.y1)/2;  i = c; j = 3;
                        }
*/
                    }
                g_tCollision = { i, j };
                SDL_SetRenderDrawColor(pSdlRenderer, 0, 0, 0, SDL_ALPHA_OPAQUE);
                RenderCircle( pSdlRenderer, 20, xd, yd );
                RenderCircle( pSdlRenderer, 21, xd, yd );
                RenderCircle( pSdlRenderer, 22, xd, yd );
                }

                SDL_RenderPresent(pSdlRenderer);

                SDL_Event event;
                while (SDL_PollEvent(&event))
                    {
                    if (event.type == SDL_QUIT)
                        {
                        done = SDL_TRUE;
                        }
                    }
                }
            }

        if (pSdlRenderer)
            {
            SDL_DestroyRenderer(pSdlRenderer);
            }
        if (pSdlWindow)
            {
            SDL_DestroyWindow(pSdlWindow);
            }
        }
    SDL_Quit();
    return 0;
    }

