#include <stack>
#include "boxes.h"
#include "ex_parse.h"
#include "notebook.h"
#include "box_convert.h"

void monobox::print(size_t indent, double offx, double offy)
{
    for (size_t i = 0; i < indent; i++)
        printf("    ");
    printf("(%f,%f) mono: size (%f,%f:%f) cursor (%d,%d)\n", offx, offy, sizex, sizey, centery, row_cursor, col_cursor);

    size_t r = 0;

    for (auto y = array.begin(); y != array.end(); ++y)
    {
        for (size_t i = 0; i <= indent; i++)
            printf("    ");

        std::cout << "line " << r << ": ";

        for (auto x = y->begin(); x != y->end(); ++x)
            std::cout << char(*x);

        std::cout << std::endl;
    }
}

visitRet monobox::visit(visitArg m)
{
    return visitret_OK;
}

boxbase * monobox::copy()
{
    return new monobox(row_cursor, col_cursor, row_cursor2, col_cursor2, array);
}

bool monobox::sort_cursor(slong & startx, slong & stopx, slong & starty, slong & stopy)
{
    starty = row_cursor;
    stopy = row_cursor2;
    startx = col_cursor;
    stopx = col_cursor2;

    if (starty < stopy)
    {
        return true;
    }
    else if (starty > stopy)
    {
        std::swap(startx, stopx);
        std::swap(starty, stopy);
        return true;
    }
    else if (startx > stopx)
    {
        std::swap(startx, stopx);
        std::swap(starty, stopy);
        return true;
    }
    else if (col_cursor < col_cursor2)
    {
        return true;
    }
    else
    {
        return false;
    }
}



void monobox::key_copy(boxbase*&b)
{
    assert(b == nullptr);

    slong startx, stopx, starty, stopy;
    if (!sort_cursor(startx, stopx, starty, stopy))
        return;

    rowbox* r = new rowbox(0, 0, 0);
    for (slong y = starty; y <= stopy; y++)
    {
        if (y > starty)
        {
            r->child.push_back(iboxarrayelem(new nullbox()));
            r->child.back().offx = 0;
            r->child.back().offy = 0;
        }

        for (slong x = (y == starty ? startx : 0); x < (y == stopy ? stopx : array[y].size()); ++x)
        {
            char16_t c = array[y][x]&65535;
            r->child.push_back(iboxarrayelem(iboximm_make(c)));
            r->child.back().offx = 0;
            r->child.back().offy = 0;
        }
    }
    b = r;
}

void monobox::key_paste(boxbase*&b)
{
    boxbase* d = nullptr;
    if (b->get_type() == BNTYPE_MONO)
    {
        monobox* B = dynamic_cast<monobox*>(b);
        for (auto y = B->array.begin(); y != B->array.end(); ++y)
        {
            for (auto x = y->begin(); x != y->end(); ++x)
                insert_char(*x);

            insert(d, insertarg_Newline);
        }
    }
    else if (b->get_type() == BNTYPE_ROW)
    {
        rowbox* B = dynamic_cast<rowbox*>(b);
        for (auto i = B->child.begin(); i != B->child.end(); ++i)
        {
            int32_t t = ibox_type(i->cibox);
            if (t >= 0)
            {
                insert_char(t);
            }
            else if (t == BNTYPE_NULLER)
            {
                insert(d, insertarg_Newline);                
            }
        }
    }
}

void monobox::delete_selection()
{
    slong startx, stopx, starty, stopy;
    if (!sort_cursor(startx, stopx, starty, stopy))
        return;
    if (stopy > starty)
    {
        array[starty].erase(array[starty].begin() + startx, array[starty].end());
        array[starty].insert(array[starty].end(), array[stopy].begin() + stopx, array[stopy].end());
        array.erase(array.begin() + starty + 1, array.begin() + stopy + 1);
    }
    else
    {
        array[starty].erase(array[starty].begin() + startx, array[starty].begin() + stopx);
    }
    col_cursor = col_cursor2 = startx;
    row_cursor = row_cursor2 = starty;
    return;
}

void monobox::insert_char(int32_t c)
{
    delete_selection();

    assert(row_cursor < array.size());
    assert(col_cursor <= array[row_cursor].size());
    
    array[row_cursor].insert(array[row_cursor].begin() + col_cursor, c);

    col_cursor++;
    col_cursor2 = col_cursor;
    row_cursor2 = row_cursor;

    int32_t * our_line = array[row_cursor].data();
    size_t our_size = array[row_cursor].size();
    int32_t r, l, cc;
    char t[40];

    // look for \[ to the left
    for (l = col_cursor - 1; l >= 0; l--)
    {
        if ((our_line[l]&65535) == '\\')
            break;
    }
    if (l < 0)
        goto alias_scan;
    if (!(l+1 < our_size && our_line[l+1]&65535 == '['))
        goto alias_scan;
    if (l > 0 && our_line[l-1]&65535 == '[')
        goto alias_scan;

    // look for ] to the right
    for (r = col_cursor - 1; r < our_size; r++)
    {
        if ((our_line[r]&65535) == '\\')
            break;
    }
    if (r >= our_size)
        goto alias_scan;

    if (r-l > 30)
        goto alias_scan;

    for (int32_t i=l; i<=r; i++)
        t[i-l] = our_line[i]&65535;
    t[r+1-l] = 0;
    cc = escapedname_to_char(t) && 65535;
    if (cc == 0)
        goto alias_scan;

    array[row_cursor].erase(array[row_cursor].begin() + l+1, array[row_cursor].begin() + r+1);
    array[row_cursor][l] = cc;
    col_cursor2 = col_cursor = l + 1;
    return;

alias_scan:

    // look for \[AliasDelimiter] to the right
    for (r = col_cursor - 1; r < our_size; r++)
    {
        if ((our_line[r]&65535) == CHAR_AliasDelimiter)
            break;
    }
    if (r >= our_size)
        return;

    // look for \[AliasDelimiter] to the left
    for (l = std::min(int32_t(col_cursor) - 1, r-1); l >= 0; l--)
    {
        if ((our_line[l]&65535) == CHAR_AliasDelimiter)
            break;
    }
    if (l<0)
        return;

    if (r-l > 30)
        return;

    for (int32_t i=l+1; i<r; i++)
        t[i-(l+1)] = our_line[i]&65535;
    t[r-(l+1)] = 0;
    cc = escape_seq_to_action(t);
    if (cc <= 0)
        return;

    array[row_cursor].erase(array[row_cursor].begin() + l+1, array[row_cursor].begin() + r+1);
    array[row_cursor][l] = cc;
    col_cursor2 = col_cursor = l + 1;
    return;
}

moveRet monobox::move(boxbase*&b, moveArg m)
{
    moveRet r = moveret_OK;
    switch (m)
    {
        case movearg_Left:
        case movearg_ShiftLeft:
        {
            if (col_cursor > 0)
            {
                col_cursor--;
            }
            else if (row_cursor > 0)
            {
                row_cursor--;
                col_cursor = array[row_cursor].size();
            }
            else
            {
                r = moveret_End;
            }

            if (m == movearg_Left)
            {
                col_cursor2 = col_cursor;
                row_cursor2 = row_cursor;                
            }
            return r;
        }

        case movearg_Right:
        case movearg_ShiftRight:
        {
            if (col_cursor < array[row_cursor].size())
            {
                col_cursor++;
            }
            else if (row_cursor + 1 < array.size())
            {
                row_cursor++;
                col_cursor = 0;
            }
            else
            {
                r = moveret_End;
            }

            if (m == movearg_Right)
            {
                col_cursor2 = col_cursor;
                row_cursor2 = row_cursor;                
            }
            return r;
        }

        case movearg_Up:
        case movearg_ShiftUp:
        {
            if (row_cursor > 0)
            {
                row_cursor--;
                col_cursor = std::min(col_cursor, (slong)array[row_cursor].size());
            }
            else
            {
                r = moveret_End;
            }

            if (m == movearg_Up)
            {
                col_cursor2 = col_cursor;
                row_cursor2 = row_cursor;
            }
            return r;
        }

        case movearg_Down:
        case movearg_ShiftDown:
        {
            if (row_cursor + 1 < array.size())
            {
                row_cursor++;
                col_cursor = std::min(col_cursor, (slong)array[row_cursor].size());
            }
            else
            {
                r = moveret_End;
            }

            if (m == movearg_Down)
            {
                col_cursor2 = col_cursor;
                row_cursor2 = row_cursor;
            }
            return r;
        }

        case movearg_Last:
        {
            row_cursor2 = row_cursor = array.size() - 1;
            col_cursor2 = col_cursor = array[row_cursor].size();
            return r;
        }
        case movearg_First:
        {
            row_cursor2 = row_cursor = 0;
            col_cursor2 = col_cursor = 0;
            return r;
        }
        case movearg_Home:
        {
            row_cursor2 = row_cursor;
            col_cursor2 = col_cursor = 0;
            return r;
        }
        case movearg_End:
        {
            row_cursor2 = row_cursor;
            col_cursor2 = col_cursor = array[row_cursor].size();
            return r;
        }
        case movearg_CtrlSpace:
        {
            if (row_cursor == row_cursor2 && col_cursor == col_cursor2)
                return moveret_End;
            else
                return moveret_OK;
        }
        case movearg_Tab:
        {
            do {
                insert_char(' ');
            } while ((col_cursor % 4) != 0);
            return r;
        }
        default:
        {
            return r;
        }
    }
}

insertRet monobox::insert(boxbase*&b, insertArg m)
{
    insertRet r = insertret_Done;
    switch (m)
    {
        case insertarg_Newline:
        {
            delete_selection();
            array.insert(array.begin() + row_cursor + 1, std::vector<int32_t>());
            for (slong x = col_cursor; x < array[row_cursor].size(); x++)
                array[row_cursor + 1].push_back(array[row_cursor][x]);
            array[row_cursor].resize(col_cursor);
            row_cursor2 = row_cursor = row_cursor + 1;
            col_cursor2 = col_cursor = 0;
            return r;
        }
        default:
        {
            return r;
        }
    }
}

removeRet monobox::remove(boxbase*&b, removeArg m)
{
    assert(b == nullptr);
    removeRet r = removeret_OK;

    if (row_cursor == row_cursor2 && col_cursor == col_cursor2)
    {
        switch (m)
        {
            case removearg_Left:
            {
                if (col_cursor > 0)
                {
                    col_cursor--;
                    array[row_cursor].erase(array[row_cursor].begin() + col_cursor);
                }
                else if (row_cursor > 0)
                {
                    col_cursor = array[row_cursor-1].size();
                    for (slong x = 0; x < array[row_cursor].size(); x++)
                        array[row_cursor - 1].push_back(array[row_cursor][x]);
                    array.erase(array.begin() + row_cursor);
                    row_cursor--;
                }
                else
                {
                    r = removeret_End;
                }

                row_cursor2 = row_cursor;
                col_cursor2 = col_cursor;
                return r;
            }
            case removearg_Right:
            {
                if (col_cursor < array[row_cursor].size())
                {
                    array[row_cursor].erase(array[row_cursor].begin() + col_cursor);
                }
                else if (row_cursor + 1 < array.size())
                {
                    for (slong x = 0; x < array[row_cursor + 1].size(); x++)
                        array[row_cursor].push_back(array[row_cursor + 1][x]);
                    array.erase(array.begin() + row_cursor + 1);
                }
                else
                {
                    r = removeret_End;
                }

                row_cursor2 = row_cursor;
                col_cursor2 = col_cursor;
                return r;
            }
            default:
            {
                assert(false);
                return r;
            }
        }
    }
    else
    {
        delete_selection();
        return removeret_OK;
    }
}

ex monobox::get_ex()
{
    std::string s;
    uex r; r.init_push_backr(gs.sym_sList.get(), array.size());
    for (auto y = array.begin(); y != array.end(); ++y)
    {
        s.clear();
        for (auto x = y->begin(); x != y->end(); ++x)
            stdstring_pushback_char16(s, char16_t(*x));
        r.push_back(emake_str(s));
    }
    return emake_node(gs.sym_sMonoBox.copy(), r.release());
}


void monobox::get_cursor(aftransform * T)
{
    T->orig_x = margin.left + col_cursor*char_sizex;
    T->orig_y = margin.top + row_cursor*line_sizey + char_centery;
    T->theta = 0;
    T->cos_theta = 1.0;
    T->sin_theta = 0.0;
}

void monobox::measure(boxmeasurearg ma)
{
    fcolor = 0;
    fs = fontsize_set_family(FONTFAMILY_LUCIDA_REG, ma.style->get_font());
    char_sizex = fontsize_char_sizex(fs, 'm');
    char_sizey = fontsize_char_sizey(fs, 'm');
    char_centery = fontsize_char_centery(fs, 'm');
    line_sizey = (17.0/16)*char_sizey;
    margin = {0.0, 0.0, 0.0, 0.0};
//    margin.scale(ma.style->nb->magnification);

    size_t m = 1;
    size_t lsize = 0;
    for (auto y = array.begin(); y != array.end(); ++y)
    {
        lsize += y->size() + 1;
        m = std::max(m, y->size());
    }
    sizex = (margin.left + margin.right) + m*char_sizex;
    sizey = (margin.top + margin.bottom + char_sizey) + (array.size() - 1)*line_sizey;
    centery = margin.top + char_centery;

    blexer L(lsize);
    for (auto y = array.begin(); y != array.end(); ++y)
    {
        for (auto x = y->begin(); x != y->end(); ++x)
        {
            L.add_char(*x);
        }
        L.add_newline();
    }

    size_t loff = 0;
    for (auto y = array.begin(); y != array.end(); ++y)
    {
        for (auto x = y->begin(); x != y->end(); ++x)
        {
            *x = (L.type[loff] << 16) + (*x & 65535);
            loff++;
        }
        loff++;
    }
}

void monobox::draw_pre(boxdrawarg da)
{
    drawtrect(da.globx, da.globx + sizex,
              da.globy, da.globy + sizey, 0x00F8F8F8, da.T);

    slong startx, stopx, starty, stopy;
    if (sort_cursor(startx, stopx, starty, stopy))
    {
        for (slong y = stopy; y > starty; y--)
        {
            if (y == stopy)
            {
                drawtrect(da.globx + margin.left + 0*char_sizex, da.globx + margin.left + stopx*char_sizex,
                          da.globy + margin.top + y*line_sizey, da.globy + margin.top + y*line_sizey + char_sizey, da.nb->cSelectionBackground, da.T);
            }
            else
            {
                drawtrect(da.globx + margin.left + 0*char_sizex, da.globx + margin.left + array[y].size()*char_sizex,
                          da.globy + margin.top + y*line_sizey, da.globy + margin.top + y*line_sizey + char_sizey, da.nb->cSelectionBackground, da.T);
            }
        }

        if (stopy == starty)
        {
            drawtrect(da.globx + margin.left + startx*char_sizex, da.globx + margin.left + stopx*char_sizex,
                      da.globy + margin.top + starty*line_sizey, da.globy + margin.top + starty*line_sizey + char_sizey, da.nb->cSelectionBackground, da.T);
        }
        else
        {
            drawtrect(da.globx + margin.left + startx*char_sizex, da.globx + margin.left + array[starty].size()*char_sizex,
                      da.globy + margin.top + starty*line_sizey, da.globy + margin.top + starty*line_sizey + char_sizey, da.nb->cSelectionBackground, da.T);
        }
    }

draw_cursor: 

    double usx = da.globx + margin.left + col_cursor*char_sizex;
    double usy = da.globy + margin.top + row_cursor*line_sizey;
    drawtline(usx, usy,
              usx, usy + char_sizey,
              0.25*fontsize_size(fs), da.nb->cCursor, da.T);
}

void monobox::draw_main(boxdrawarg da)
{
    if (da.dflags & DFLAG_SELECTION)
    {
        for (slong y = 0; y < array.size(); y++)
        {
            for (slong x = 0; x < array[y].size(); x++)
            {
                drawtchar(fs, array[y][x]&65535,
                          char_sizex, char_sizey, da.globx + margin.left + x*char_sizex, da.globy + margin.top + y*line_sizey,
                          da.nb->cSelectionForeground, da.T);
            }
        }
        return;
    }

    slong startx, stopx, starty, stopy;

    if (da.dflags & DFLAG_IGNORESEL)
    {
        startx = stopx = starty = stopy = 0;
    }
    else
    {
        sort_cursor(startx, stopx, starty, stopy);
    }

    if (da.dflags & DFLAG_SCOLOR)
    {
        for (slong y = 0; y < array.size(); y++)
        {
            for (slong x = 0; x < array[y].size(); x++)
            {
                int32_t child_type = array[y][x];
                uint32_t color = da.nb->pallet1[(child_type>>16)&255];
                if (((child_type>>16)&255) == lextype_symbol_1st ||
                    ((child_type>>16)&255) == lextype_pattern_1st ||
                    ((child_type>>16)&255) == lextype_symbol)
                {
                    color = da.nb->pallet1[lextype_MAX + ((child_type>>24)&255)];
                }

                if (y == starty)
                {
                    if (x >= startx)
                    {
                        if (y != stopy || x < stopx)
                            color = da.nb->cSelectionForeground;
                    }
                }
                else if (y == stopy)
                {
                    if (x < stopx)
                        color = da.nb->cSelectionForeground;
                }
                else if (starty < y && y < stopy)
                {
                    color = da.nb->cSelectionForeground;
                }

                drawtchar(fs, child_type&65535,
                          char_sizex, char_sizey, da.globx + margin.left + x*char_sizex, da.globy + margin.top + y*line_sizey,
                          color, da.T);
            }
        }
    }
    else
    {
        for (slong y = 0; y < array.size(); y++)
        {
            for (slong x = 0; x < array[y].size(); x++)
            {
                uint32_t color = fcolor;

                if (y == starty)
                {
                    if (x >= startx)
                    {
                        if (y != stopy || x < stopx)
                            color = da.nb->cSelectionForeground;
                    }
                }
                else if (y == stopy)
                {
                    if (x < stopx)
                        color = da.nb->cSelectionForeground;
                }
                else if (starty < y && y < stopy)
                {
                    color = da.nb->cSelectionForeground;
                }

                drawtchar(fs, array[y][x]&65535,
                          char_sizex, char_sizey, da.globx + margin.left + x*char_sizex, da.globy + margin.top + y*line_sizey,
                          fcolor, da.T);
            }
        }
    }
}

void monobox::draw_post(boxdrawarg da)
{

}

