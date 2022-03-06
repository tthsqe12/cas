#include <stack>
#include "boxes.h"
#include "ex_parse.h"
#include "notebook.h"
#include "box_convert.h"

void gridbox::print(size_t indent, double offx, double offy)
{
    for (size_t i = 0; i < indent; i++)
        printf("    ");
    printf("(%f,%f) grid: size (%f,%f:%f) cursor (%d,%d)\n", offx, offy, sizex, sizey, centery, col_cursor, row_cursor);

    for (auto i = array.begin(); i != array.end(); ++i)
        for (auto j = i->begin(); j != i->end(); ++j)
            j->cbox->print(indent + 1, j->offx, j->offy);
}

visitRet gridbox::visit(visitArg m)
{
    for (auto i = array.begin(); i != array.end(); ++i)
        for (auto j = i->begin(); j != i->end(); ++j)
            j->cbox->visit(m);
    return visitret_OK;
}

boxbase * gridbox::copy()
{
    gridbox * r = new gridbox(array, row_cursor, col_cursor);
    for (auto i = r->array.begin(); i != r->array.end(); ++i)
        for (auto j = i->begin(); j != i->end(); ++j)
            j->cbox = dynamic_cast<rowbox*>(j->cbox->copy());
    return r;
}

void gridbox::key_copy(boxbase*&b)
{
    if (row_cursor2 == row_cursor && col_cursor2 == col_cursor)
    {
        return;
    }
    else
    {
        slong xa = std::min(col_cursor, col_cursor2);
        slong xb = std::max(col_cursor, col_cursor2);
        slong ya = std::min(row_cursor, row_cursor2);
        slong yb = std::max(row_cursor, row_cursor2);

        std::vector<std::vector<rowboxarrayelem>> newarray;
        for (slong i = ya; i <= yb; i++)
        {
            std::vector<rowboxarrayelem> newrow;
            for (slong j = xa; j <= xb; j++)
            {
                newrow.push_back(rowboxarrayelem(dynamic_cast<rowbox*>(array[i][j].cbox->copy())));
            }
            newarray.push_back(newrow);
        }

        assert(b == nullptr);
        b = new gridbox(newarray, 0, 0);
    }
}

void gridbox::key_paste(boxbase*&b)
{
    if (row_cursor2 == row_cursor && col_cursor2 == col_cursor)
    {
        array[row_cursor][col_cursor].cbox->key_paste(b);
    }
    return;
}

void gridbox::insert_char(int32_t c)
{
    array[row_cursor][col_cursor].cbox->insert_char(c);
}

moveRet gridbox::move(boxbase*&b, moveArg m)
{
    assert(b == nullptr);
    moveRet r;
    rowbox * us = array[row_cursor][col_cursor].cbox;

    switch (m)
    {
        case movearg_Left:
        {
            row_cursor2 = row_cursor;
            col_cursor2 = col_cursor;
            if (us->is_selected_placeholder())
                r = moveret_End;
            else
                r = us->move(b, m);
            assert(r == moveret_OK || r == moveret_End);
            if (r == moveret_End)
            {
                if (col_cursor > 0)
                {
                    col_cursor2 = col_cursor = col_cursor - 1;
                    array[row_cursor][col_cursor].cbox->move(b, movearg_Last);
                    return moveret_OK;
                }
                else if (row_cursor > 0)
                {
                    row_cursor2 = row_cursor = row_cursor - 1;
                    col_cursor2 = col_cursor = array[row_cursor].size() - 1;
                    array[row_cursor][col_cursor].cbox->move(b, movearg_Last);
                    return moveret_OK;
                }
            }
            return r;
        }
        case movearg_ShiftLeft:
        {
            if (row_cursor2 == row_cursor && col_cursor2 == col_cursor)
            {
                r = us->move(b, m);
                assert(r == moveret_OK || r == moveret_End);
                if (r == moveret_End && col_cursor > 0)
                {
                    col_cursor -= 1;
                    return moveret_OK;
                }
                else
                {
                    return r;
                }
            }
            else if (col_cursor > 0)
            {
                col_cursor -= 1;
                return moveret_OK;
            }
            else
            {
                return moveret_End;
            }
        }
        case movearg_Right:
        {
            row_cursor2 = row_cursor;
            col_cursor2 = col_cursor;
            if (us->is_selected_placeholder())
                r = moveret_End;
            else
                r = us->move(b, m);
            assert(r == moveret_OK || r == moveret_End);
            if (r == moveret_End)
            {
                if (col_cursor < array[row_cursor].size() - 1)
                {
                    col_cursor2 = col_cursor = col_cursor + 1;
                    array[row_cursor][col_cursor].cbox->move(b, movearg_First);
                    return moveret_OK;
                }
                else if (row_cursor < array.size() - 1)
                {
                    row_cursor2 = row_cursor = row_cursor + 1;
                    col_cursor2 = col_cursor = 0;
                    array[row_cursor][col_cursor].cbox->move(b, movearg_First);
                    return moveret_OK;
                }
            }
            return r;
        }
        case movearg_ShiftRight:
        {
            if (row_cursor2 == row_cursor && col_cursor2 == col_cursor)
            {
                r = us->move(b, m);
                assert(r == moveret_OK || r == moveret_End);
                if (r == moveret_End && col_cursor + 1 < array[row_cursor].size())
                {
                    col_cursor += 1;
                    return moveret_OK;
                }
                else
                {
                    return r;
                }
            }
            else if (col_cursor + 1 < array[row_cursor].size())
            {
                col_cursor += 1;
                return moveret_OK;
            }
            else
            {
                return moveret_End;
            }
        }
        case movearg_Last:
        {
            row_cursor = row_cursor2 = array.size() - 1;
            col_cursor = col_cursor2 = array[row_cursor].size() - 1;
            array[row_cursor][col_cursor].cbox->move(b, movearg_Last);
            return moveret_OK;
        }
        case movearg_First:
        {
            row_cursor = row_cursor2 = 0;
            col_cursor = col_cursor2 = 0;
            array[row_cursor][col_cursor].cbox->move(b, movearg_First);
            return moveret_OK;
        }
        case movearg_Home:
        case movearg_End:
        {
            row_cursor2 = row_cursor;
            col_cursor2 = col_cursor;
            r = us->move(b, m);
            assert(r == moveret_OK || r == moveret_End);
            return r;
        }
        case movearg_Switch:
        {
            if (row_cursor2 == row_cursor && col_cursor2 == col_cursor)
            {
                r = array[row_cursor][col_cursor].cbox->move(b, m);
                assert(r == moveret_OK || r == moveret_End);
                return r;
            }
            else
            {
                return moveret_End;
            }
        }
        case movearg_Up:
        {
            r = us->move(b, m);
            assert(r == moveret_OK || r == moveret_End);
            if (r != moveret_OK && row_cursor > 0)
            {
                row_cursor--;
                col_cursor = std::min(col_cursor, (slong)array[row_cursor].size() - 1);
                r = moveret_OK;
            }
            row_cursor2 = row_cursor;
            col_cursor2 = col_cursor;
            return r;
        }
        case movearg_Down:
        {
            r = us->move(b, m);
            assert(r == moveret_OK || r == moveret_End);
            if (r != moveret_OK && row_cursor + 1 < array.size())
            {
                row_cursor += 1;
                col_cursor = std::min(col_cursor, (slong)array[row_cursor].size() - 1);
                r = moveret_OK;
            }
            row_cursor2 = row_cursor;
            col_cursor2 = col_cursor;
            return r;    
        }
        case movearg_ShiftDown:
        {
            if (row_cursor2 == row_cursor && col_cursor2 == col_cursor)
            {
                r = us->move(b, m);
                assert(r == moveret_OK || r == moveret_End);
                if (r == moveret_End && row_cursor + 1 < array.size())
                {
                    row_cursor += 1;
                    col_cursor = std::min(col_cursor, (slong)array[row_cursor].size() - 1);
                    return moveret_OK;
                }
                else
                {
                    return r;
                }
            }
            else if (row_cursor + 1 < array.size())
            {
                row_cursor += 1;
                col_cursor = std::min(col_cursor, (slong)array[row_cursor].size() - 1);
                return moveret_OK;
            }
            else
            {
                return moveret_End;
            }
        }
        case movearg_ShiftUp:
        {
            if (row_cursor2 == row_cursor && col_cursor2 == col_cursor)
            {
                r = us->move(b, m);
                assert(r == moveret_OK || r == moveret_End);
                if (r == moveret_End && row_cursor > 0)
                {
                    row_cursor -= 1;
                    col_cursor = std::min(col_cursor, (slong)array[row_cursor].size() - 1);
                    return moveret_OK;
                }
                else
                {
                    return r;
                }
            }
            else if (row_cursor > 0)
            {
                row_cursor -= 1;
                col_cursor = std::min(col_cursor, (slong)array[row_cursor].size() - 1);
                return moveret_OK;
            }
            else
            {
                return moveret_End;
            }
        }

        case movearg_CtrlSpace:
        {
            if (row_cursor == row_cursor2 && col_cursor == col_cursor2)
            {
                r = us->move(b, m);
                assert(r == moveret_OK || r == moveret_End);
                return r;
            }
            else
            {
                return moveret_OK;
            }
        }
        default:
        {
            assert(false);
            return moveret_OK;
        }
    }
}

rowbox * new_placeholder()
{
    rowbox* newstuff = new rowbox(2, 0,1);
    newstuff->child[0].cibox = iboximm_make(CHAR_Placeholder);
    newstuff->child[1].cibox.ptr = new nullbox();
    return newstuff;
}

insertRet gridbox::insert(boxbase*&b, insertArg m)
{
    rowbox * us = array[row_cursor][col_cursor].cbox;

    if (m == insertarg_GridRow && us->has_cursor())
    {
        size_t n = array[row_cursor].size();
        row_cursor++;
        array.insert(array.begin() + row_cursor, std::vector<rowboxarrayelem>());
        array[row_cursor].push_back(rowboxarrayelem(new_placeholder()));
        for (size_t i = 1; i < n; i++)
            array[row_cursor].push_back(rowboxarrayelem(new_placeholder()));
        col_cursor = 0;
        row_cursor2 = row_cursor;
        col_cursor2 = col_cursor;
        return insertret_Done;
    }
    else if (m == insertarg_GridCol && us->has_cursor())
    {
        col_cursor++;
        for (auto & i : array)
        {
            while (i.size() < col_cursor)
                i.push_back(rowboxarrayelem(new_placeholder()));
            i.insert(i.begin() + col_cursor, rowboxarrayelem(new_placeholder()));
        }
        row_cursor = 0;
        row_cursor2 = row_cursor;
        col_cursor2 = col_cursor;
        return insertret_Done;
    }
    else
    {
        return us->insert(b, m);
    }
}

removeRet gridbox::remove(boxbase*&b, removeArg m)
{
    assert(b == nullptr);
    removeRet r;

    if (row_cursor == row_cursor2 && col_cursor == col_cursor2)
    {
        rowbox * us = array[row_cursor][col_cursor].cbox;
        if (us->is_selected_placeholder())
        {
            if (m == removearg_Left)
            {
                if (col_cursor > 0)
                {
                    col_cursor2 = col_cursor = col_cursor - 1;
                    array[row_cursor][col_cursor].cbox->move(b, movearg_Last);
                    return removeret_OK;
                }
                for (auto j = array[row_cursor].begin(); j != array[row_cursor].end(); ++j)
                {
                    if (j->cbox->is_empty() || j->cbox->is_placeholder())
                        continue;

                    if (row_cursor > 0)
                    {
                        row_cursor2 = row_cursor = row_cursor - 1;
                        col_cursor2 = col_cursor = array[row_cursor].size() - 1;
                        array[row_cursor][col_cursor].cbox->move(b, movearg_Last);
                        return removeret_OK;
                    }
                    else
                    {
                        return removeret_End;
                    }
                }
                array.erase(array.begin() + row_cursor);
                if (array.empty())
                {
                    b = nullptr;
                    return removeret_Replace;
                }
                else
                {
                    if (row_cursor > 0)
                        row_cursor--;
                    row_cursor2 = row_cursor;
                    col_cursor2 = col_cursor = array[row_cursor].size() - 1;
                    array[row_cursor][col_cursor].cbox->move(b, movearg_Last);
                    return removeret_OK;
                }
            }
            else
            {
                if (col_cursor + 1 < array[row_cursor].size())
                {
                    col_cursor2 = col_cursor = col_cursor + 1;
                    array[row_cursor][col_cursor].cbox->move(b, movearg_First);
                    return removeret_OK;
                }
                for (auto j = array[row_cursor].begin(); j != array[row_cursor].end(); ++j)
                {
                    if (j->cbox->is_empty() || j->cbox->is_placeholder())
                        continue;

                    if (row_cursor + 1 < array[row_cursor].size())
                    {
                        row_cursor2 = row_cursor = row_cursor + 1;
                        col_cursor2 = col_cursor = 0;
                        array[row_cursor][col_cursor].cbox->move(b, movearg_First);
                        return removeret_OK;
                    }
                    else
                    {
                        return removeret_End;
                    }
                }
                array.erase(array.begin() + row_cursor);
                if (array.empty())
                {
                    b = nullptr;
                    return removeret_Replace;
                }
                else
                {
                    if (row_cursor > 0)
                        row_cursor--;
                    row_cursor2 = row_cursor;
                    col_cursor2 = col_cursor = 0;
                    array[row_cursor][col_cursor].cbox->move(b, movearg_First);
                    return removeret_OK;
                }
            }
        }
        else
        {
            r = us->remove(b, m);
            assert(r == removeret_End || r == removeret_OK);
            us->made_into_placeholder();
            return removeret_OK;
        }
    }
    else
    {
        slong xa = std::min(col_cursor, col_cursor2);
        slong xb = std::max(col_cursor, col_cursor2);
        slong ya = std::min(row_cursor, row_cursor2);
        slong yb = std::max(row_cursor, row_cursor2);

        for (slong i = ya; i <= yb; i++)
        {
            for (slong j = xa; j <= xb; j++)
            {
                if (j < array[i].size())
                    array[i][j].cbox->set_placeholder();
            }
        }

        for (slong i = yb; i >= ya; i--)
        {
            bool ok = true;
            for (auto j = array[i].begin(); j != array[i].end(); ++j)
            {
                if (j->cbox->is_empty() || j->cbox->is_placeholder())
                    continue;

                ok = false;
                break;
            }
            if (ok)
                array.erase(array.begin() + i);
        }

        for (slong j = xb; j >= xa; j--)
        {
            bool ok = true;
            for (auto i = array.begin(); i != array.end(); ++i)
            {
                if (j >= i->size() || i->at(j).cbox->is_empty() || i->at(j).cbox->is_placeholder())
                    continue;

                ok = false;
                break;
            }
            if (ok)
            {
                for (auto i = array.begin(); i != array.end(); ++i)
                {
                    if (j >= i->size())
                        continue;

                    i->erase(i->begin() + j);
                }
            }
        }

        if (array.empty())
        {
            b = nullptr;
            return removeret_Replace;
        }
        else
        {
            row_cursor = std::min(ya, (slong)array.size() - 1);
            col_cursor = std::min(xa, (slong)array[row_cursor].size() - 1);
            col_cursor2 = col_cursor;
            row_cursor2 = row_cursor;
            return removeret_OK;
        }
    }
}

ex gridbox::get_ex()
{
    std::vector<uex> l(1);
    l.back().init_push_backr(gs.sym_sList.get(), array.size());
    for (auto i = array.begin(); i != array.end(); ++i)
    {
        uex r; r.init_push_backr(gs.sym_sList.get(), i->size());
        for (auto j = i->begin(); j != i->end(); ++j)
            r.push_back(j->cbox->get_ex());
        l.back().push_back(r.release());
    }

    if (!eis_sym(alignment.get(), gs.sym_sNull.get()))
    {
        l.push_back(emake_node(gs.sym_sRule.copy(), gs.sym_sGridBoxAlignment.copy(),
                                                    alignment.copy()));
    }

    return emake_node(gs.sym_sGridBox.copy(), l);
}



void gridbox::get_cursor(aftransform * T)
{
    rowboxarrayelem e = array[row_cursor][col_cursor];
    e.cbox->get_cursor(T);
    T->orig_x += e.offx;
    T->orig_y += e.offy;
}

void gridbox::measure(boxmeasurearg ma)
{
    uint32_t fs = ma.style->get_font();

    ulong row_count = array.size();
    ulong col_count = 0;

    std::vector<double> max_width, acc_width;
    std::vector<double> max_above(row_count, 0);
    std::vector<double> max_below(row_count, 0);

    ma.style->new_options();
    ma.style->push_opt(opt_FontSize, emake_double(29.0/32*ma.style->get_fontsize()));
    ma.style->init_const1();
    ma.level++;

    for (ulong j = 0; j < array.size(); j++)
    {
        col_count = std::max(col_count, array[j].size());
        for (ulong i = 0; i < array[j].size(); i++)
        {
            array[j][i].cbox->measure(ma);
            while (max_width.size() <= i)
                max_width.push_back(0);
            max_width[i] = std::max(max_width[i], array[j][i].cbox->sizex);
            max_above[j] = std::max(max_above[j], array[j][i].cbox->centery);
            max_below[j] = std::max(max_below[j], array[j][i].cbox->sizey - array[j][i].cbox->centery);
        }
    }

    ma.style->pop_options();

    /* read alignment options */
    std::vector<int32_t> row_option(row_count, 1);
    std::vector<int32_t> col_option(col_count, 1);
    er t = alignment.get();
    if (ehas_head_sym(t, gs.sym_sList.get()))
    {
        for (ulong i = 0; i < elength(t); i++)
        {
            er ti = echild(t,i+1);
            if (ehas_head_sym_length(ti, gs.sym_sRule.get(), 2))
            {
                er o = echild(ti,2);
                if (0 == estr_cmp(echild(ti,1), gs.strColumns.get()))
                {
                    if (ehas_head_sym(o, gs.sym_sList.get()))
                    {
                        ulong h = 0;
                        for (ulong j = 0; j < elength(o); j++)
                        {
                            if (h >= col_count)
                                break;

                            if (ehas_head_sym(echild(o,j+1), gs.sym_sList.get()))
                            {
                                for (ulong k = 0; ; k++)
                                {
                                    if (h >= col_count)
                                        break;

                                    if (k >= elength(echild(o,j+1)))
                                        k -= elength(echild(o,j+1));
                                    col_option[h++] = (eis_sym(echild(o,j+1,k+1), gs.sym_sLeft.get())) ? 0 :
                                                      (eis_sym(echild(o,j+1,k+1), gs.sym_sRight.get())) ? 2 : 1;
                                }
                            }
                            else
                            {
                                col_option[h++] = (eis_sym(echild(o,j+1), gs.sym_sLeft.get())) ? 0 :
                                                  (eis_sym(echild(o,j+1), gs.sym_sRight.get())) ? 2 : 1;
                            }
                        }
                    }
                }
                else if (eis_str(echild(ti,1), "Rows"))
                {
                    
                }
            }
        }
    }

    double f = fontsize_size(fs);
    double accx = f*GRID_EXTRAX1;
    for (ulong i = 0; i < max_width.size(); i++)
    {
        acc_width.push_back(accx);
        accx += max_width[i] + f*GRID_EXTRAX2;
    }
    accx += f*(GRID_EXTRAX3 - GRID_EXTRAX2);

    double accy = f*GRID_EXTRAY1;
    for (ulong j = 0; j < array.size(); j++)
    {
        for (ulong i = 0; i < array[j].size(); i++)
        {
            array[j][i].offx = acc_width[i] + col_option[i]*(max_width[i] - array[j][i].cbox->sizex)/2;
            array[j][i].offy = accy + (max_above[j] - array[j][i].cbox->centery);
        }
        accy += max_above[j] + max_below[j] + f*GRID_EXTRAY2;
    }
    accy += f*(GRID_EXTRAY3 - GRID_EXTRAY2);

    sizex = accx;
    sizey = accy;
    centery = sizey/2;
}



void gridbox::draw_pre(boxdrawarg da)
{
    if (row_cursor2 == row_cursor && col_cursor2 == col_cursor)
    {
        array[row_cursor][col_cursor].cbox->draw_pre(boxdrawarg(da, array[row_cursor][col_cursor].offx,
                                                                    array[row_cursor][col_cursor].offy));
    }
    else
    {
        slong xa = std::min(col_cursor, col_cursor2);
        slong xb = std::max(col_cursor, col_cursor2);
        slong ya = std::min(row_cursor, row_cursor2);
        slong yb = std::max(row_cursor, row_cursor2);

        for (slong i = ya; i <= yb; i++)
        {
            for (slong j = xa; j <= xb; j++)
            {
                double child_sizex = array[i][j].cbox->sizex;
                double child_sizey = array[i][j].cbox->sizey;
                double usx = da.globx + array[i][j].offx;
                double usy = da.globy + array[i][j].offy;
                drawtrect(usx, usx + child_sizex,
                          usy, usy + child_sizey, da.nb->cSelectionBackground, da.T);
            }
        }
    }
}

void gridbox::draw_main(boxdrawarg da)
{
    if (da.dflags & DFLAG_SELECTION)
    {
        for (ulong i = 0; i < array.size(); i++)
        {
            for (ulong j = 0; j < array[i].size(); j++)
            {
                array[i][j].cbox->draw_main(boxdrawarg(da, array[i][j].offx, array[i][j].offy, DFLAG_SELECTION));
                
            }
        }

    }
    else if (da.dflags & DFLAG_IGNORESEL)
    {
        for (ulong i = 0; i < array.size(); i++)
        {
            for (ulong j = 0; j < array[i].size(); j++)
            {
                array[i][j].cbox->draw_main(boxdrawarg(da, array[i][j].offx, array[i][j].offy, 0));
                
            }
        }

    }
    else if (row_cursor2 == row_cursor && col_cursor2 == col_cursor)
    {
        for (ulong i = 0; i < array.size(); i++)
        {
            for (ulong j = 0; j < array[i].size(); j++)
            {
                array[i][j].cbox->draw_main(boxdrawarg(da, array[i][j].offx, array[i][j].offy,
                                         (i == row_cursor && j == col_cursor) ? 0 : DFLAG_IGNORESEL));
                
            }
        }
    }    
    else
    {
        slong xa = std::min(col_cursor, col_cursor2);
        slong xb = std::max(col_cursor, col_cursor2);
        slong ya = std::min(row_cursor, row_cursor2);
        slong yb = std::max(row_cursor, row_cursor2);
        for (ulong i = 0; i < array.size(); i++)
        {
            for (ulong j = 0; j < array[i].size(); j++)
            {
                array[i][j].cbox->draw_main(boxdrawarg(da, array[i][j].offx, array[i][j].offy,
                                         (ya <= i && i <= yb && xa <= j && j <= xb) ? DFLAG_SELECTION : DFLAG_IGNORESEL));
                
            }
        }
    }
}

void gridbox::draw_post(boxdrawarg da)
{
    if (row_cursor2 == row_cursor && col_cursor2 == col_cursor)
    {
        array[row_cursor][col_cursor].cbox->draw_post(boxdrawarg(da,
                                                        array[row_cursor][col_cursor].offx,
                                                        array[row_cursor][col_cursor].offy));
    }
}
