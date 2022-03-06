#include "timing.h"
#include "font.h"
#include "graphics.h"
#include "notebook.h"
#include "box_lex.h"
#include "graphics.h"
#include "svg.h"
#include "eval.h"
#include "ex_print.h"

void cpcolorp(double*a, double *b) {a[0]=b[0]; a[1]=b[1]; a[2]=b[2];}



void drawchar(uint32_t fs, char16_t c, int sizey, int offx, int offy, double * color)
{
    uint32_t size = fs&65535;
    uint32_t family = fs >> 16;
//printf("drawing char: family = %d size = %d   offx = %d offy = %d\n", family, size, offx, offy);
/*
    if (family == FONTFAMILY_SERIFREG)
    {
        glf_serifreg.drawchar(&glb_image, glb_image.serifreg_todraw, c, size, sizey, offx, offy, color);
    }
    else if (family == FONTFAMILY_SERIFBOLD)
    {
        glf_serifbold.drawchar(&glb_image, glb_image.serifbold_todraw, c, size, sizey, offx, offy, color);
    }
    else if (family == FONTFAMILY_SANSREG)
    {
        glf_sansreg.drawchar(&glb_image, glb_image.sansreg_todraw, c, size, sizey, offx, offy, color);
    }
    else if (family == FONTFAMILY_SANSBOLD)
    {
        glf_sansbold.drawchar(&glb_image, glb_image.sansbold_todraw, c, size, sizey, offx, offy, color);
    }
*/
}

int32_t chardrawcount;

void drawtchar(uint32_t fs, char16_t c, double sizex, double sizey, double offx, double offy, uint32_t color, aftransform * T)
{
    int32_t size = fontsize_size(fs);
    uint32_t family = fontsize_family(fs);
    int32_t ci = glb_fonts[family].get_characteridx(c);
    affTran S(T->cos_theta, -T->sin_theta,
              T->sin_theta, T->cos_theta,
              T->orig_x + T->cos_theta*offx - T->sin_theta*offy,
              T->orig_y + T->sin_theta*offx + T->cos_theta*offy);
    if (ci >= 0)
    {
        renderChar(&glb_image, color, glb_chars.real_data[ci], S, sizex, sizey);
    }
    else
    {
        renderRect(&glb_image, color, S, 0, sizex, 0, sizey);
    }
}


void drawbtchar(double blend, uint32_t fs, char16_t c, double sizex, double sizey, double offx, double offy, uint32_t color, aftransform * T)
{
    int32_t size = fontsize_size(fs);
    uint32_t family = fontsize_family(fs);
    int32_t ci = glb_fonts[family].get_characteridx(c);
    affTran S(T->cos_theta, -T->sin_theta,
              T->sin_theta, T->cos_theta,
              T->orig_x + T->cos_theta*offx - T->sin_theta*offy,
              T->orig_y + T->sin_theta*offx + T->cos_theta*offy);
    if (ci >= 0)
    {
        renderBlendChar(&glb_image, blend, color, glb_chars.real_data[ci], S, sizex, sizey);
    }
    else
    {
        renderBlendRect(&glb_image, blend, color, S, 0, sizex, 0, sizey);
    }
}


void drawtline(double X1, double Y1, double X2, double Y2, double e, uint32_t color, aftransform*T)
{
    renderLine(&glb_image,
                T->orig_x + T->cos_theta*(X1) - T->sin_theta*(Y1),
                T->orig_y + T->sin_theta*(X1) + T->cos_theta*(Y1),
                T->orig_x + T->cos_theta*(X2) - T->sin_theta*(Y2),
                T->orig_y + T->sin_theta*(X2) + T->cos_theta*(Y2),
            e, color);
}


void drawtlines(double * coords, size_t nlines, double e, uint32_t color, aftransform*T)
{
    for (size_t i = 0; i < nlines; i++)
    {
        drawtline(coords[2*i+0], coords[2*i+1], coords[2*i+2], coords[2*i+3], e, color, T);
    }
}

void drawtrect(double minx, double maxx, double miny, double maxy, uint32_t color, aftransform*T)
{
    affTran S(T->cos_theta, -T->sin_theta,
              T->sin_theta, T->cos_theta,
              T->orig_x,
              T->orig_y);
    renderRect(&glb_image, color, S, minx, maxx, miny, maxy);
}

void drawbtrect(double blend, double minx, double maxx, double miny, double maxy, uint32_t color, aftransform*T)
{
    affTran S(T->cos_theta, -T->sin_theta,
              T->sin_theta, T->cos_theta,
              T->orig_x,
              T->orig_y);
    renderBlendRect(&glb_image, blend, color, S, minx, maxx, miny, maxy);
}



#if 0
void notebook::_draw_cursor_selection(int globx, int globy)
{
    sel.clear();

    aftransform T;
    get_cursor_parentT(T, globx, globy);

    if (selection.type == SELTYPE_NONE)
    {
    }
    else if (selection.type == SELTYPE_ROW)
    {
    }
    else if (selection.type == SELTYPE_GRID)
    {
    }
    else if (selection.type == SELTYPE_CELL)
    {
        assert(!selection.data.empty());
        std::vector<int> tstack;
        tstack.resize(_depth());
        for (size_t i = 0; i < _depth(); i++)
        {
            tstack[i] = cursor1[i].idx;
        }
        _push_cell_selection(bfrom_node(root), 0, tstack, 0, selection.data);
        // dont draw cursor if non trivial selection
        for (size_t i = 0; i < _depth(); i++)
        {
            assert(i < selection.data.size());
            if (cursor1[i].idx != selection.data[i])
            {
                return;
            }
        }
    }
    else
    {
        assert(false);
    }

    /* draw cursor now */

    box us = _us();
    boxnode * parent = bto_node(_parent());
    int32_t pi = _pi();

    if (boxnode_type(parent) == BNTYPE_ROOT || boxnode_type(parent) == BNTYPE_CELLGROUP)
    {
        drawtline(0,                parent->array[pi].offy,
                  glb_image.sizex,  parent->array[pi].offy,
                  zoomint*0.03, cCursor, &T);
    }
    else if (boxnode_type(parent) == BNTYPE_STUB)
    {
/*
        drawtline(0,                parent->array[pi].offy,
                  glb_image.sizex,  parent->array[pi].offy,
                  zoomint*0.03, cCursor, &T);
*/
    }
    else
    {
        assert(boxnode_type(parent) == BNTYPE_ROW);

        {
            int32_t idx  = pi;
            std::vector<int32_t> ostack;
            std::vector<int32_t> commapos;
            int32_t match = 0;
            int32_t right_pi = -1;
            int32_t right_x = 0;
            ostack.clear();
            for (int32_t count = 5000; count >= 0; count--)
            {
                int32_t t = btype(parent->array[idx].child);
                if (t >= 0)
                {
                    int32_t y, x = t & 0x0FFFFFF;
                    if (x == ',' + 65536*lextype_comma && ostack.empty())
                    {
                        commapos.push_back(idx);
                    }
                    else if (   (y =  '(' + 65536*lextype_parenth_open,
                                 x == ')' + 65536*lextype_parenth_close)
                             || (y =  '{' + 65536*lextype_curly_open,
                                 x == '}' + 65536*lextype_curly_close)
                             || (y =  '[' + 65536*lextype_bracket_open,
                                 x == ']' + 65536*lextype_bracket_close)
                             || (y =  CHAR_LeftDoubleBracket + 65536*lextype_bracket_open,
                                 x == CHAR_RightDoubleBracket + 65536*lextype_bracket_close))
                    {
                        if (ostack.empty())
                        {
                            right_pi =idx;
                            right_x = x;
                            break;
                        }
                        if (y == ostack.back())
                        {
                            ostack.pop_back();
                        }
                    }
                    else if (   x == '(' + 65536*lextype_parenth_open
                             || x == '{' + 65536*lextype_curly_open
                             || x == '[' + 65536*lextype_bracket_open
                             || x == CHAR_LeftDoubleBracket + 65536*lextype_bracket_open)
                    {
                        ostack.push_back(x);
                    }
                }
                ++idx;
                if (idx >= boxnode_len(parent))
                {
                    break;
                }
            }

            idx = pi;
            match = 0;
            int32_t left_pi = -1;
            int32_t left_x = 0;
            ostack.clear();
            for (int32_t count = 5000; count >= 0; count--)
            {
                --idx;
                if (idx < 0)
                {
                    break;
                }
                int32_t t = btype(parent->array[idx].child);
                if (t >= 0)
                {
                    int32_t y, x = t & 0x0FFFFFF;
                    if (x == ',' + 65536*lextype_comma && ostack.empty())
                    {
                        commapos.push_back(idx);
                    }
                    else if (   (y =  ')' + 65536*lextype_parenth_close,
                                 x == '(' + 65536*lextype_parenth_open)
                             || (y =  '}' + 65536*lextype_curly_close,
                                 x == '{' + 65536*lextype_curly_open)
                             || (y =  ']' + 65536*lextype_bracket_close,
                                 x == '[' + 65536*lextype_bracket_open)
                             || (y =  CHAR_RightDoubleBracket + 65536*lextype_bracket_close,
                                 x == CHAR_LeftDoubleBracket + 65536*lextype_bracket_open))
                    {
                        if (ostack.size() < 1)
                        {
                            left_pi = idx;
                            left_x = x;
                            break;
                        }
                        if (y == ostack.back())
                        {
                            ostack.pop_back();
                        }
                    }
                    else if (   x == ')' + 65536*lextype_parenth_close
                             || x == '}' + 65536*lextype_curly_close
                             || x == ']' + 65536*lextype_bracket_close
                             || x == CHAR_RightDoubleBracket + 65536*lextype_bracket_close)
                    {
                        ostack.push_back(x);
                    }
                }
            }
/*
            rasterfont * f = fontint_to_fontp(bnode_extra0(parent));
            int32_t default_sizey = f->chars[99].sizey;
            int32_t default_centery = f->chars[99].centery;
*/
            box sib;
            int32_t sib_offx, sib_offy;
            int32_t sib_type, sib_sizex, sib_sizey, sib_centery;

            uint32_t bracketcolor = cBracketMismatch;
            if (   (    left_x == '[' + 65536*lextype_bracket_open
                    && right_x == ']' + 65536*lextype_bracket_close)
                || (    left_x == '(' + 65536*lextype_parenth_open
                    && right_x == ')' + 65536*lextype_parenth_close)
                || (    left_x == '{' + 65536*lextype_curly_open
                    && right_x == '}' + 65536*lextype_curly_close)
                || (    left_x == CHAR_LeftDoubleBracket + 65536*lextype_bracket_open
                    && right_x == CHAR_RightDoubleBracket + 65536*lextype_bracket_close))
            {
                bracketcolor = cBracketMatch;
                for (size_t i = 0; i < commapos.size(); i++)
                {
                    sib = parent->array[commapos[i]].child;
                    sib_offx = parent->array[commapos[i]].offx;
                    sib_offy = parent->array[commapos[i]].offy;
                    bget_header(sib_type, sib_sizex, sib_sizey, sib_centery, sib);
                    drawtrect(sib_offx, sib_offx + sib_sizex, sib_offy + sib_sizey/2, sib_offy + sib_sizey, bracketcolor, &T);
                }
            }
            if (left_x != 0)
            {
                sib = parent->array[left_pi].child;
                sib_offx = parent->array[left_pi].offx;
                sib_offy = parent->array[left_pi].offy;
                bget_header(sib_type, sib_sizex, sib_sizey, sib_centery, sib);
//std::cout << "highlight left " << std::endl;
//                bitmap_draw_rectangle(image, 0x00c0ff90, parentx + sib_offx, parenty + sib_offy, sib_sizex, sib_sizey);
//                bitmap_draw_rchar(f->chars[sib_type&255], image, parentx + sib_offx, parenty + sib_offy, sib_sizey, bracketcolor);

                drawtrect(sib_offx, sib_offx + sib_sizex, sib_offy, sib_offy + sib_sizey, bracketcolor, &T);

            }
            if (right_x != 0)
            {
                sib = parent->array[right_pi].child;
                sib_offx = parent->array[right_pi].offx;
                sib_offy = parent->array[right_pi].offy;
                bget_header(sib_type, sib_sizex, sib_sizey, sib_centery, sib);
//std::cout << "highlight right " << std::endl;
//                bitmap_draw_rectangle(image, 0x00c0ff90, parentx + sib_offx, parenty + sib_offy, sib_sizex, sib_sizey);
//                bitmap_draw_rchar(f->chars[sib_type&255], image, parentx + sib_offx, parenty + sib_offy, sib_sizey, bracketcolor);
                drawtrect(sib_offx, sib_offx + sib_sizex, sib_offy, sib_offy + sib_sizey, bracketcolor, &T);
            }
        }



        uint32_t fs = boxnode_extra0(parent);
        int32_t default_sizey = fontsize_default_sizey(fs);
        int32_t default_centery = fontsize_default_centery(fs);
        int32_t us_type, us_sizex, us_sizey, us_centery;
        bget_header(us_type, us_sizex, us_sizey, us_centery, us);
        int32_t cursor_y;
        int32_t cursor_sizey;
        if (pi > 0 && btype(boxnode_child(parent, pi - 1)) != BNTYPE_NULLER)
        {
            int32_t prev_type, prev_sizex, prev_sizey, prev_centery;
            int32_t us_offy = boxnode_offy(parent, pi);
            int32_t prev_offy = boxnode_offy(parent, pi - 1);
            bget_header(prev_type, prev_sizex, prev_sizey, prev_centery, boxnode_child(parent, pi - 1));
            if (us_type != BNTYPE_NULLER)
            {
                cursor_y = (us_centery - prev_centery)/2;
                cursor_sizey = (prev_sizey + us_sizey)/2;
            }
            else
            {
                cursor_y = us_centery - prev_centery;
                cursor_sizey = prev_sizey;
            }
        }
        else
        {
            if (us_type != BNTYPE_NULLER)
            {
                cursor_y = 0;
                cursor_sizey = us_sizey;
            }
            else
            {
                cursor_y = us_centery - default_centery;
                cursor_sizey = default_sizey;
            }
        }

        if (us_type == BNTYPE_ROT)
        {
            int32_t usx = parent->array[pi].offx;
            int32_t usy = parent->array[pi].offy;
            double Ax = usx;
            double Ay = usy;
            double Bx = usx + us_sizex;
            double By = usy;
            double Cx = usx;
            double Cy = usy + us_sizey;
            double Dx = usx + us_sizex;
            double Dy = usy + us_sizey;
            drawtline(Ax, Ay, Bx, By, 3.0 + fontsize_size(fs)*0.02, glb_image.fcolor, &T);
            drawtline(Bx, By, Dx, Dy, 3.0 + fontsize_size(fs)*0.02, glb_image.fcolor, &T);
            drawtline(Dx, Dy, Cx, Cy, 3.0 + fontsize_size(fs)*0.02, glb_image.fcolor, &T);
            drawtline(Cx, Cy, Ax, Ay, 4.0 + fontsize_size(fs)*0.04, cCursor, &T);
        }
        else
        {
            drawtline(parent->array[pi].offx, parent->array[pi].offy + cursor_y,
                      parent->array[pi].offx, parent->array[pi].offy + cursor_y + cursor_sizey,
                      6.0 + fontsize_size(fs)*0.04, cCursor, &T);
/*
            int32_t pk = pi;
            int32_t max_y = 0;
            if (!(bis_char(parent->array[pk].child)
                  && isletterchar(char16_t(bchar_type(parent->array[pk].child)))))
            {
                while (pk > 0 && bis_char(parent->array[pk-1].child)
                              && isletterchar(char16_t(bchar_type(parent->array[pk-1].child))))
                {
                    int32_t child_type, child_sizex, child_sizey, child_centery;
                    bget_header(child_type, child_sizex, child_sizey, child_centery, parent->array[pk-1].child);
                    pk--;
                    max_y = std::max(max_y, parent->array[pk].offy + child_sizey);
                }
            }

            if (pk < pi - 1)
            {
                std::string s;
                for (int32_t k = pk; k < pi; k++)
                {
                    stdstring_pushback_char16(s, bchar_type(parent->array[k].child)&65535);
                }
std::cout << "prefix " << s << std::endl;                    

                // Prefix search
                auto set_range = gs.char_set.equal_prefix_range(s);

                drawtline(parent->array[pk].offx, max_y,
                          parent->array[pi].offx, max_y,
                          6.0 + fontsize_size(fs)*0.06, cCursor, &T);

                int32_t ox, oy = max_y;

                for(auto it = set_range.first; it != set_range.second; ++it)
                {
                    ox = parent->array[pk].offx;
                    it.key(s);
                    size_t sn = s.size();
                    size_t si;
                    int32_t max_above = 0, max_below = 0;
                    si = 0;
                    while (si < sn)
                    {
                        char16_t c;
                        si += readonechar16(c, reinterpret_cast<const unsigned char *>(s.data()) + si);
                        int32_t c_fsizey = fontsize_char_sizey(fs, c);
                        int32_t c_fcentery = fontsize_char_centery(fs, c);
                        max_above = std::max(max_above, c_fcentery);
                        max_below = std::max(max_below, c_fsizey - c_fcentery);
                    }
                    si = 0;
                    while (si < sn)
                    {
                        char16_t c;
                        si += readonechar16(c, reinterpret_cast<const unsigned char *>(s.data()) + si);
                        int32_t c_fsizex = fontsize_char_sizex(fs, c);
                        int32_t c_fsizey = fontsize_char_sizey(fs, c);
                        int32_t c_fcentery = fontsize_char_centery(fs, c);

                        drawtchar(fs, c, c_fsizex, c_fsizey, ox, oy + max_above - c_fcentery, 0x00808080, &T);
                        ox += c_fsizex;
                    }
                    oy += max_above + max_below + 1;
                }
            }
*/
        }
    }
}








void notebook::_draw_post(int globx, int globy)
{
    sel.clear();

    aftransform T;
    get_cursor_parentT(T, globx, globy);

    if (selection.type != SELTYPE_NONE)
    {
        return;
    }

    box us = _us();
    boxnode * parent = bto_node(_parent());
    int32_t pi = _pi();

    if (boxnode_type(parent) == BNTYPE_ROOT || boxnode_type(parent) == BNTYPE_CELLGROUP)
    {
        return;
    }
    else if (boxnode_type(parent) == BNTYPE_STUB)
    {
        return;
    }
    else
    {
        assert(boxnode_type(parent) == BNTYPE_ROW);

        uint32_t fs = boxnode_extra0(parent);
        int32_t default_sizey = fontsize_default_sizey(fs);
        int32_t default_centery = fontsize_default_centery(fs);
        int32_t us_type, us_sizex, us_sizey, us_centery;
        bget_header(us_type, us_sizex, us_sizey, us_centery, us);
        int32_t cursor_y;
        int32_t cursor_sizey;

        {
            int32_t pk = pi;
            int32_t max_y = 0;
            if (!(bis_char(parent->array[pk].child)
                  && isletterchar(char16_t(bchar_type(parent->array[pk].child)))))
            {
                while (pk > 0 && bis_char(parent->array[pk-1].child)
                              && isletterchar(char16_t(bchar_type(parent->array[pk-1].child))))
                {
                    int32_t child_type, child_sizex, child_sizey, child_centery;
                    bget_header(child_type, child_sizex, child_sizey, child_centery, parent->array[pk-1].child);
                    pk--;
                    max_y = std::max(max_y, parent->array[pk].offy + child_sizey);
                }
            }

            if (pk < pi - 1)
            {
                std::string s;
                for (int32_t k = pk; k < pi; k++)
                {
                    stdstring_pushback_char16(s, bchar_type(parent->array[k].child)&65535);
                }
//std::cout << "prefix " << s << std::endl;                    

                // Prefix search
                auto set_range = gs.char_set.equal_prefix_range(s);

                int32_t oox, ox, oy = max_y;

                for(auto it = set_range.first; it != set_range.second; ++it)
                {
                    oox = parent->array[pk].offx;
                    it.key(s);
                    size_t sn = s.size();
                    size_t si;
                    int32_t max_above = 0, max_below = 0;
                    si = 0;
                    ox = oox;
                    while (si < sn)
                    {
                        char16_t c;
                        si += readonechar16(c, reinterpret_cast<const unsigned char *>(s.data()) + si);
                        int32_t c_fsizex = fontsize_char_sizex(fs, c);
                        int32_t c_fsizey = fontsize_char_sizey(fs, c);
                        int32_t c_fcentery = fontsize_char_centery(fs, c);
                        max_above = std::max(max_above, c_fcentery);
                        max_below = std::max(max_below, c_fsizey - c_fcentery);
                        ox += c_fsizex;
                    }

                    drawbtrect(0.75, oox, ox, oy, oy + max_above + max_below, 0x00f0f0f0, &T);

                    si = 0;
                    ox = oox;
                    while (si < sn)
                    {
                        char16_t c;
                        si += readonechar16(c, reinterpret_cast<const unsigned char *>(s.data()) + si);
                        int32_t c_fsizex = fontsize_char_sizex(fs, c);
                        int32_t c_fsizey = fontsize_char_sizey(fs, c);
                        int32_t c_fcentery = fontsize_char_centery(fs, c);

                        drawbtchar(0.75, fs, c, c_fsizex, c_fsizey, ox, oy + max_above - c_fcentery, 0x00101010, &T);
                        ox += c_fsizex;
                    }
                    oy += max_above + max_below;
                }
            }
        }
    }
}







bool notebook::_handle_cellclick(std::vector<int>&s, box Us, int level, int X, int Y)
{
//printf("cellclick(%d)(%d,%d)\n",level,X,Y);// boxnode_print(NULL,Us,0);

    int cellbracket_w = glb_dingbat.get_char_sizex(DINGBAT_CELLGEN, zoomint);
    int cellbracket_h = glb_dingbat.get_char_sizey(DINGBAT_CELLGEN, zoomint);

    if (Y < 0 || Y >= bnode_sizey(Us))
    {
        return false;
    }
    else if (level > 0 && X >= glb_image.sizex - cellbracket_w*(2*level+1)/2
                       && X <  glb_image.sizex - cellbracket_w*(2*level-1)/2)
    {
//printf("found\n");// boxnode_print(NULL,Us,0);

            boxnode * us = bto_node(Us);

            assert(us->header.type == BNTYPE_CELL || us->header.type == BNTYPE_CELLGROUP);

            // find starting corrodinate y and height for bracket
            int32_t y     = 0;
            int32_t sizey = us->header.sizey;
            int32_t firstoffy = 0;
            int32_t lastoffy = 0;
            boxnode * first = us;
            while (first->header.type == BNTYPE_CELLGROUP)
            {
                firstoffy += first->array[0].offy;
                first = bto_node(first->array[0].child);
            }
            boxnode *  last = us;
            while (last->header.type == BNTYPE_CELLGROUP)
            {
                if (last->extra1 & BNFLAG_OPEN)
                {
                    lastoffy +=  last->array[last->len-1].offy;
                    last = bto_node(last->array[last->len-1].child);
                }
                else
                {
                    lastoffy +=  last->array[0].offy;
                    last = bto_node(last->array[0].child);
                }
            }
            firstoffy += first->array[0].offy; first = bto_node(first->array[0].child);
            int32_t firstsizey = first->header.sizey;
            lastoffy += last->array[0].offy; last = bto_node(last->array[0].child);
            int32_t lastsizey = last->header.sizey;

            int delta = 3;
            firstoffy -= delta; firstsizey += 2*delta;
            lastoffy -= delta; lastsizey += 2*delta;
            // make sure first/last sizey is big enough for dingbat
            if (firstsizey < cellbracket_h)
            {
                delta = cellbracket_h - firstsizey;
                firstoffy -= delta/2; firstsizey += delta;
            }
            if (lastsizey < cellbracket_h)
            {
                delta = cellbracket_h - lastsizey;
                lastoffy -= delta/2; lastsizey += delta;
            }

            // we have the location and size of the bracket now
            y = 0 + firstoffy;
            sizey = lastsizey + lastoffy - firstoffy;

        return y<=Y && Y<y+sizey;
    }
    else if (bnode_type(Us) == BNTYPE_CELL)
    {
        return false;
    }
    else
    {
        assert(bnode_type(Us) == BNTYPE_CELLGROUP);
        for (int32_t i = (bnode_extra1(Us) & BNFLAG_OPEN) ? bnode_len(Us) - 1 : 0; i >= 0; i--)
        {
            if (_handle_cellclick(s, bnode_child(Us,i), level + 1, X, Y - bnode_offy(Us,i)))
            {
                s.push_back(i);
                return true;
            }
        }
        return false;
    }
}

bool notebook::handle_cellclick(int X, int Y, bool doubleclicked)
{
    std::vector<int> s;
    for (int i = root->len - 2; i >= 0; i--)
    {
        s.clear();
        if (_handle_cellclick(s, root->array[i].child, 1, X, Y - root->array[i].offy))
        {
            cursor1.clear();
            cursor1.push_back(cursorentry(root));
//printf("a down %d\n",i);
            _down1(i);
            for (size_t j = s.size(); j > 0; j--) {_down1(s[j-1]);}
            if (doubleclicked && btype(_us()) == BNTYPE_CELLGROUP)
            {
                bto_node(_us())->extra1 ^= BNFLAG_OPEN;
                bto_node(_us())->extra1 &= ~BNFLAG_MEASURED;
            }
            while (btype(_us()) == BNTYPE_CELLGROUP) {_down1(0);}

//printf("first:\n"); print();

            selection.type = SELTYPE_CELL;
            selection.data.resize(_depth());
            for (size_t j = 0; j < _depth(); j++) {selection.data[j] = cursor1[j].idx;}

            cursor1.clear();
            cursor1.push_back(cursorentry(root));
            _down1(i);
//printf("b down %d\n",i);
            for (size_t j = s.size(); j > 0; j--) {_down1(s[j-1]);}
            while (_pi() + 1 >= bnode_len(_parent())) {_up1();}
            _right1();
            while (btype(_us()) == BNTYPE_CELLGROUP) {_down1(0);}

//printf("last:\n"); print();


            return true;
        }
    }
    return false;
}


void notebook::bitmap_draw_cellgroup_bracket(boxnode * us, int level, double globx, double globy, bool selected)
{
//std::cout << "bitmap_draw_cellgroup_selection golbx = " << globx << " globy = " << globy << std::endl;
//boxnode_print(bfrom_ptr(NULL), bfrom_node(us), 0);

    assert(us->header.type == BNTYPE_CELL || us->header.type == BNTYPE_CELLGROUP);

    /* find starting corrodinate y and height for bracket */
    double y     = globy;
    double sizey = us->header.sizey;
    double firstoffy = 0;
    double lastoffy = 0;
    boxnode * first = us;


    while (first->header.type == BNTYPE_CELLGROUP)
    {
        firstoffy += first->array[0].offy;
        first = bto_node(first->array[0].child);
    }
    boxnode *  last = us;
    while (last->header.type == BNTYPE_CELLGROUP)
    {
        if (last->extra1 & BNFLAG_OPEN)
        {
            lastoffy +=  last->array[last->len-1].offy;
            last = bto_node(last->array[last->len-1].child);
        }
        else
        {
            lastoffy +=  last->array[0].offy;
            last = bto_node(last->array[0].child);
        }
    }
    firstoffy += first->array[0].offy; first = bto_node(first->array[0].child);
    double firstsizey = first->header.sizey;
    lastoffy += last->array[0].offy; last = bto_node(last->array[0].child);
    double lastsizey = last->header.sizey;

std::cout << "------------" << std::endl;
std::cout << "firstoffy: " << firstoffy << std::endl;
std::cout << " lastoffy: " << lastoffy << std::endl;
std::cout << "firstsizey: " << firstsizey << std::endl;
std::cout << " lastsizey: " << lastsizey << std::endl;


    double cellbracket_w = glb_dingbat.get_char_sizex(DINGBAT_CELLGEN, zoomint);
    double cellbracket_h = glb_dingbat.get_char_sizey(DINGBAT_CELLGEN, zoomint);

    double delta = 3;
    firstoffy -= delta; firstsizey += 2*delta;
    lastoffy -= delta; lastsizey += 2*delta;
    /* make sure first/last sizey is big enough for dingbat */
    if (firstsizey < cellbracket_h)
    {
        delta = cellbracket_h - firstsizey;
        firstoffy -= delta/2; firstsizey += delta;
    }
    if (lastsizey < cellbracket_h)
    {
        delta = cellbracket_h - lastsizey;
        lastoffy -= delta/2; lastsizey += delta;
    }

    /* we have the location and size of the bracket now */
    y = globy + firstoffy;
    sizey = lastsizey + lastoffy - firstoffy;

    int x = glb_image.sizex - cellbracket_w*(2*level+1)/2;
    uint32_t color = cCellBracket;
    if (selected)
    {        
        color = cSelectionForeground;

        aftransform T;
        T.orig_x = 0.0;
        T.orig_y = 0.0;
        T.theta = 0;
        T.rix = 2.0/glb_image.sizex;
        T.riy = 2.0/glb_image.sizey;
        T.cos_theta = 1.0;
        T.sin_theta = 0.0;
        drawtrect(x, x + cellbracket_w, y, y + sizey, cSelectionBackground, &T);
    }

    char16_t c;
    if (us->header.type == BNTYPE_CELLGROUP)
    {
        c = (us->extra1 & BNFLAG_OPEN) ? DINGBAT_CELLGEN : DINGBAT_CELLCLOSED;
    }
    else
    {
        c = us->extra0 == CELLTYPE_INPUT   ? DINGBAT_CELLINPUT :
            us->extra0 == CELLTYPE_PRINT   ? DINGBAT_CELLOUTPUT :
            us->extra0 == CELLTYPE_OUTPUT  ? DINGBAT_CELLOUTPUT :
            us->extra0 == CELLTYPE_MESSAGE ? DINGBAT_CELLMESSAGE :
                                             DINGBAT_CELLTEXT;
    }
    glb_dingbat.draw_char(&glb_image, c, color, affTran(1,0,0,1,x,y), cellbracket_w, sizey);
}


void notebook::_push_cell_selection(box Us, size_t i, const std::vector<int>&a, size_t j, const std::vector<int>&b)
{
/*
std::cout << "_draw_cell_selection i = " << i << " j = " << j << std::endl;
boxnode_print(bfrom_ptr(NULL), bfrom_node(us), 0);
*/
    if (j < i)
    {
        boxnode * us = bto_node(Us);
        assert(j == 0);
        // draw everything after and including a
        if (i == a.size())
        {
            sel.push_back(us);
        }
        else
        {
            if (a[i] == 0)
            {
                sel.push_back(us);
            }
            else
            {
                for (int32_t k = a[i] + 1; k < us->len; k++)
                {
                    sel.push_back(bto_node(us->array[k].child));
                }
                _push_cell_selection(us->array[a[i]].child, i+1,a,0,b);
                return;
            }
        }
    }
    else if (i < j)
    {
        assert(i == 0);
        // draw everything strictly before b
        if (j < b.size())
        {
            boxnode * us = bto_node(Us);
            for (int32_t k = 0; k < b[j]; k++)
            {
                sel.push_back(bto_node(us->array[k].child));
            }
            _push_cell_selection(us->array[b[j]].child, 0,a,j+1,b);
            return;
        }
    }
    else
    {
        assert(i == j);
        if (i >= a.size() || j >= b.size())
        {
            // push nothing
            return;
        }
        boxnode * us = bto_node(Us);
        if (a[i] > b[j])
        {
            for (int32_t k = b[j] + 1; k < a[i]; k++)
            {
                sel.push_back(bto_node(us->array[k].child));
            }
            _push_cell_selection(us->array[b[j]].child, j+1,b,0,a);
            _push_cell_selection(us->array[a[i]].child, 0,b,i+1,a);
        }
        else if (a[i] < b[j])
        {
            for (int32_t k = a[i] + 1; k < b[j]; k++)
            {
                sel.push_back(bto_node(us->array[k].child));
            }
            _push_cell_selection(us->array[a[i]].child, i+1,a,0,b);
            _push_cell_selection(us->array[b[j]].child, 0,a,j+1,b);
        }
        else
        {
            _push_cell_selection(us->array[a[i]].child, i+1,a,j+1,b);
        }
    }
}
#endif


void _draw_cellgroup_bracket(boxbase * us, boxdrawarg da, bool selected)
{
//std::cout << "_draw_cellgroup_bracket called" << std::endl;

    assert(us->get_type() == BNTYPE_CELL || us->get_type() == BNTYPE_CELLGROUP);

    /* find starting corrodinate y and height for bracket */
    double y     = da.globy;
    double sizey = us->sizey;
    double firstoffy = 0;
    double lastoffy = 0;
    boxbase * first = us;

    while (first->get_type() == BNTYPE_CELLGROUP)
    {
        firstoffy += offyat(first, 0);
        first = childat(first, 0);
    }
    boxbase* last = us;

    while (last->get_type() == BNTYPE_CELLGROUP)
    {
/*
        if (last->extra1 & BNFLAG_OPEN)
        {
*/
            lastoffy += offyat(last, childlen(last) - 1);
            last = childat(last, childlen(last) - 1);
/*
        }
        else
        {
//std::cout << " reading from first !!!!!!!!!!!!!!!!!" << std::endl;
            lastoffy +=  last->array[0].offy;
            last = bto_node(last->array[0].child);
        }
*/
    }

    cellbox* ffirst = dynamic_cast<cellbox*>(first);
    cellbox* llast = dynamic_cast<cellbox*>(last);

    
    firstoffy += ffirst->bracket_offy;
    double firstsizey = ffirst->bracket_sizey;
    lastoffy += llast->bracket_offy;
    double lastsizey = llast->bracket_sizey;

    double cellbracket_w = glb_dingbat.get_char_sizex(DINGBAT_CELLGEN, da.nb->magnification*11.0);
    double cellbracket_h = glb_dingbat.get_char_sizey(DINGBAT_CELLGEN, da.nb->magnification*11.0);

    double delta = 3;
    firstoffy -= delta; firstsizey += 2*delta;
    lastoffy -= delta; lastsizey += 2*delta;

    /* make sure first/last sizey is big enough for dingbat */
    if (firstsizey < cellbracket_h)
    {
        delta = cellbracket_h - firstsizey;
        firstoffy -= 0.5*delta; firstsizey += delta;
    }
    if (lastsizey < cellbracket_h)
    {
        delta = cellbracket_h - lastsizey;
        lastoffy -= 0.5*delta; lastsizey += delta;
    }

    /* we have the location and size of the bracket now */
    y = da.globy + firstoffy;
    sizey = lastsizey + lastoffy - firstoffy;

    int x = glb_image.sizex - cellbracket_w*0.5*(2*da.level+1);
    uint32_t color = da.nb->cCellBracket;
    if (selected)
    {        
        color = da.nb->cSelectionForeground;
        glb_dingbat.draw_char(&glb_image, DINGBAT_CELLHIGHLIGHT, da.nb->cSelectionBackground, affTran(1,0,0,1,x,y), cellbracket_w, sizey);
    }

    char16_t c = DINGBAT_CELLGEN;
/*
    if (us->get_type() == BNTYPE_CELLGROUP)
    {
        c = (true) ? DINGBAT_CELLGEN : DINGBAT_CELLCLOSED;
    }
    else
    {
        cellType ct = dynamic_cast<cellbox*>(us)->celltype;
        c = ct == cellt_Input   ? DINGBAT_CELLINPUT :
            ct == cellt_Print   ? DINGBAT_CELLOUTPUT :
            ct == cellt_Output  ? DINGBAT_CELLOUTPUT :
            ct == cellt_Message ? DINGBAT_CELLMESSAGE :
                                  DINGBAT_CELLTEXT;
    }
*/
    glb_dingbat.draw_char(&glb_image, c, color, affTran(1,0,0,1,x,y), cellbracket_w, sizey);
}



/*
    else if (us->header.type == BNTYPE_BUTTON)
    {
        assert(us->len == 1);
        drawtrect(globx, globx + us->header.sizex, globy, globy + us->header.sizey, RGB_COLOR(220, 220, 220), T);
        for (int i = 0; i < 1; i++)
        {
            _draw_bitmap(bto_node(us->array[i].child), level+1, globx + us->array[i].offx, globy + us->array[i].offy, dflags, T);
        }
    }
    else if (us->header.type == BNTYPE_GRID)
    {
        assert(us->len == us->extra0*us->extra1);
        if (selection.type == SELTYPE_GRID && us == bto_node(_gparent()))
        {
            assert(selection.data.size() == 1);
            int gridi = _gpi();

            int a_x = selection.data[0] % us->extra0;
            int a_y = selection.data[0] / us->extra0;        
            int b_x = gridi % us->extra0;
            int b_y = gridi / us->extra0;
            if (a_x > b_x) {std::swap(a_x, b_x);}
            if (a_y > b_y) {std::swap(a_y, b_y);}

            int * minx = (int *) malloc(us->extra0*sizeof(int));
            int * maxx = (int *) malloc(us->extra0*sizeof(int));
            int * miny = (int *) malloc(us->extra1*sizeof(int));
            int * maxy = (int *) malloc(us->extra1*sizeof(int));
            for (int i = 0; i < us->extra0; i++)
            {
                minx[i] = 2000000000;
                maxx[i] = 0;
            }
            for (int j = 0; j < us->extra1; j++)
            {
                miny[j] = 2000000000;
                maxy[j] = 0;
            }
            for (int j = 0; j < us->extra1; j++) {
                for (int i = 0; i < us->extra0; i++) {
                    boxheader * child = bto_ptr(us->array[us->extra0*j + i].child);
                    minx[i] = std::min(minx[i], us->array[us->extra0*j + i].offx);
                    maxx[i] = std::max(maxx[i], us->array[us->extra0*j + i].offx + child->sizex);
                    miny[j] = std::min(miny[j], us->array[us->extra0*j + i].offy);
                    maxy[j] = std::max(maxy[j], us->array[us->extra0*j + i].offy + child->sizey);
                }
            }

            int mx = globx + (a_x > 0 ? (maxx[a_x-1] + minx[a_x])/2 : 0);
            int my = globy + (a_y > 0 ? (maxy[a_y-1] + miny[a_y])/2 : 0);
            int dx = (b_x + 1 < us->extra0 ? (maxx[b_x] + minx[b_x+1])/2 : us->header.sizex) - (a_x > 0 ? (maxx[a_x-1] + minx[a_x])/2 : 0);
            int dy = (b_y + 1 < us->extra1 ? (maxy[b_y] + miny[b_y+1])/2 : us->header.sizey) - (a_y > 0 ? (maxy[a_y-1] + miny[a_y])/2 : 0);

            drawtrect(mx, mx + dx, my, my + dy, cSelectionBackground, T);

            free(maxy);
            free(miny);
            free(maxx);
            free(minx);

            for (int j = 0; j < us->extra1; j++)
            {
                for (int i = 0; i < us->extra0; i++)
                {
                    _draw_bitmap(bto_node(us->array[us->extra0*j + i].child), level+1,
                                     globx + us->array[us->extra0*j + i].offx,
                                     globy + us->array[us->extra0*j + i].offy,
                                     dflags | ((a_x <= i && i <= b_x && a_y <= j && j <= b_y) ? DFLAG_SELECTION : 0), T);
                }
            }
        }
        else
        {
            for (int j = 0; j < us->extra1; j++)
            {
                for (int i = 0; i < us->extra0; i++)
                {
                    _draw_bitmap(bto_node(us->array[us->extra0*j + i].child), level+1,
                                    globx + us->array[us->extra0*j + i].offx,
                                    globy + us->array[us->extra0*j + i].offy,
                                    dflags, T);
                }
            }
        }
    }
    else if (us->header.type == BNTYPE_TABVIEW)
    {
        drawtrect(globx, globx + us->header.sizex, globy, globy + us->header.sizey, RGB_COLOR(240, 240, 240), T);
        for (int32_t j = 0; j < us->len; j += 2)
        {
            _draw_bitmap(bto_node(us->array[j].child), level+1,
                             globx + us->array[j].offx,
                             globy + us->array[j].offy,
                             dflags, T);
        }

        int32_t j = us->extra0;
        assert(j < us->len/2);
            _draw_bitmap(bto_node(us->array[2*j+1].child), level+1,
                             globx + us->array[2*j+1].offx,
                             globy + us->array[2*j+1].offy,
                             dflags, T);
    }
*/
