#include <string>
#include <iostream>
#include <functional> // std::optional


#include "grundpunkt.h"
#include "te.h"

using namespace std::string_literals;

VEbenenLagen g_vEbenen;     // E1,  E2,  E3  = 3 homologe Lagen einer Ebene
VPolDreieck  g_vPoldreieck; // P12, P13, P23 = 3 Polpunkte zu den o.g.Ebenenlagen
VGelenke     g_vGelenke;
SCollision   g_tCollision;

CGrundpunkt  g_oGrundpunkt ( {-1, -1} );
VGrundpunkte g_vGrundpunkte;  // A123, B123

// using TRenderItem = std::map<std::string, std::string>;
// using TRenderData = std::multimap<std::string, TRenderItem>;

template<typename T>
    auto L( T const & P1, T const & P2 )
	{
	return sqrt( pow((P1.x-P2.x), 2) + pow((P1.y-P2.y), 2) );
	}

void ExportSCAD( VEbenenLagen const & crvEL,
	         VGrundpunkte const & crvGP )
    {
    auto const A0 = crvGP[0].G0();
    auto const B0 = crvGP[1].G0();
    auto const E1 = crvEL[0];
    auto const E2 = crvEL[1];
    auto const E3 = crvEL[2];

    auto const a1 = crvGP[0].GPoint(0);
    auto const b1 = crvGP[1].GPoint(0);
    auto const a2 = crvGP[0].GPoint(1);
    auto const b2 = crvGP[1].GPoint(1);
    auto const a3 = crvGP[0].GPoint(2);
    auto const b3 = crvGP[1].GPoint(2);

    auto const GL = L(A0, B0);
    auto const AL = L(A0, a1);
    auto const BL = L(B0, b1);
    auto const CL = L(a1, b1);

    TRenderItem oSubM{};
    TRenderData oData{};

    oSubM.emplace("A0.x", std::to_string(A0.x));
    oSubM.emplace("A0.y", std::to_string(A0.y));
    oSubM.emplace("B0.x", std::to_string(B0.x));
    oSubM.emplace("B0.y", std::to_string(B0.y));

    oSubM.emplace("a1.x", std::to_string(a1.x));
    oSubM.emplace("a1.y", std::to_string(a1.y));
    oSubM.emplace("b1.x", std::to_string(b1.x));
    oSubM.emplace("b1.y", std::to_string(b1.y));

    oSubM.emplace("a2.x", std::to_string(a2.x));
    oSubM.emplace("a2.y", std::to_string(a2.y));
    oSubM.emplace("b2.x", std::to_string(b2.x));
    oSubM.emplace("b2.y", std::to_string(b2.y));

    oSubM.emplace("a3.x", std::to_string(a3.x));
    oSubM.emplace("a3.y", std::to_string(a3.y));
    oSubM.emplace("b3.x", std::to_string(b3.x));
    oSubM.emplace("b3.y", std::to_string(b3.y));

    oSubM.emplace("GL",  std::to_string(GL));
    oSubM.emplace("AL",  std::to_string(AL));
    oSubM.emplace("BL",  std::to_string(BL));
    oSubM.emplace("CL",  std::to_string(CL));

    oData.emplace( "Getriebe"s, oSubM );
    oSubM.clear();

    oSubM.emplace("E1P1.x",  std::to_string(E1.x1));
    oSubM.emplace("E1P1.y",  std::to_string(E1.y1));
    oSubM.emplace("E1P2.x",  std::to_string(E1.x2));
    oSubM.emplace("E1P2.y",  std::to_string(E1.y2));

    oSubM.emplace("E2P1.x",  std::to_string(E2.x1));
    oSubM.emplace("E2P1.y",  std::to_string(E2.y1));
    oSubM.emplace("E2P2.x",  std::to_string(E2.x2));
    oSubM.emplace("E2P2.y",  std::to_string(E2.y2));

    oSubM.emplace("E3P1.x",  std::to_string(E3.x1));
    oSubM.emplace("E3P1.y",  std::to_string(E3.y1));
    oSubM.emplace("E3P2.x",  std::to_string(E3.x2));
    oSubM.emplace("E3P2.y",  std::to_string(E3.y2));

    oData.emplace( "Ebene"s, oSubM );
    oSubM.clear();

    Cte ote(oData, "3LagenSynthese.tmpl", "../templates/");
    std::cout << ote << std::endl;
    }



/******************************************************************************************

*/


int main(int argc, char* argv[])
    {
    MakeCircle();

    if (SDL_Init(SDL_INIT_VIDEO) == 0)
        {
        SDL_Window   * pSdlWindow   = nullptr;
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
			    if ( g_vGrundpunkte.size() == 2 ) g_vGrundpunkte.erase(g_vGrundpunkte.begin());
			    if ( g_vGrundpunkte.size() < 2 )
				{
				g_vGrundpunkte.emplace_back(g_oGrundpunkt);
				if ( g_vGrundpunkte.size() == 2 )
				    {
				    ExportSCAD( g_vEbenen, g_vGrundpunkte);
				    }
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

                SDL_SetRenderDrawColor(pSdlRenderer, 0xc0, 0xc0, 0xc0, SDL_ALPHA_OPAQUE);
                for (int l = 0; l < 27; ++l)
                    {
                    SDL_RenderDrawLine(pSdlRenderer,  l*100, 0, l*100, 26*100 );
                    SDL_RenderDrawLine(pSdlRenderer,  0, l*100, 26*100, l*100 );
                    }

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

