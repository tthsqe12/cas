#include <stack>
#include "boxes.h"
#include "ex_parse.h"
#include "notebook.h"
#include "box_convert.h"

void fractionbox::print(size_t indent, double offx, double offy)
{
    for (size_t i = 0; i < indent; i++)
        printf("    ");
    printf("(%f,%f) fraction: size(%f,%f:%f) cursor(%d)\n", offx, offy, sizex, sizey, sizey, cursor);

    num.cbox->print(indent + 1, num.offx, num.offy);
    den.cbox->print(indent + 1, den.offx, den.offy);
}

visitRet fractionbox::visit(visitArg m)
{
    num.cbox->visit(m);
    den.cbox->visit(m);
    return visitret_OK;
}

boxbase * fractionbox::copy()
{
    rowbox * newnum = dynamic_cast<rowbox*>(num.cbox->copy());
    rowbox * newden = dynamic_cast<rowbox*>(den.cbox->copy());
    fractionbox * r = new fractionbox(newnum, newden, cursor);
    r->num.offx = num.offx;
    r->num.offy = num.offy;
    r->den.offx = den.offx;
    r->den.offy = den.offy;
    r->fs = fs;
    return r;    
}

void fractionbox::key_copy(boxbase*&b)
{
    if (cursor == 0)
    {
        num.cbox->key_copy(b);
    }
    else
    {
        assert(cursor == 1);
        den.cbox->key_copy(b);
    }
}

void fractionbox::key_paste(boxbase*&b)
{
    if (cursor == 0)
    {
        num.cbox->key_paste(b);
    }
    else
    {
        assert(cursor == 1);
        den.cbox->key_paste(b);
    }
}

void fractionbox::insert_char(int32_t c)
{
    switch (cursor)
    {
        case 0:
            num.cbox->insert_char(c);
            return;
        case 1:
            den.cbox->insert_char(c);
            return;
        default:
            assert(false);
    }
}

moveRet fractionbox::move(boxbase*&b, moveArg m)
{
    moveRet r;
    switch (m)
    {
        case movearg_Left:
        {
            if (cursor == 0)
            {
                r = num.cbox->move(b, movearg_Left);
                assert(r == moveret_OK || r == moveret_End);
                return r;
            }
            else
            {
                assert(cursor == 1);
                r = den.cbox->move(b, movearg_Left);
                assert(r == moveret_OK || r == moveret_End);
                if (r == moveret_End)
                {
                    r = num.cbox->move(b, movearg_Last);
                    assert(r == moveret_OK);
                    cursor = 0;
                }
                return moveret_OK;
            }
        }
        case movearg_ShiftLeft:
        {
            if (cursor == 0)
            {
                r = num.cbox->move(b, movearg_ShiftLeft);
                assert(r == moveret_OK || r == moveret_End);
                return r;
            }
            else
            {
                assert(cursor == 1);
                r = den.cbox->move(b, movearg_ShiftLeft);
                assert(r == moveret_OK || r == moveret_End);
                return r;
            }
        }
        case movearg_Right:
        {
            if (cursor == 1)
            {
                r = den.cbox->move(b, movearg_Right);
                assert(r == moveret_OK || r == moveret_End);
                return r;
            }
            else
            {
                assert(cursor == 0);
                r = num.cbox->move(b, movearg_Right);
                assert(r == moveret_OK || r == moveret_End);
                if (r == moveret_End)
                {
                    r = den.cbox->move(b, movearg_First);
                    assert(r == moveret_OK);
                    cursor = 1;
                }
                return moveret_OK;
            }
        }
        case movearg_ShiftRight:
        {
            if (cursor == 0)
            {
                r = num.cbox->move(b, movearg_ShiftRight);
                assert(r == moveret_OK || r == moveret_End);
                return r;
            }
            else
            {
                assert(cursor == 1);
                r = den.cbox->move(b, movearg_ShiftRight);
                assert(r == moveret_OK || r == moveret_End);
                return r;
            }
        }
        case movearg_Last:
        {
            r = den.cbox->move(b, movearg_Last);
            assert(r == moveret_OK);
            cursor = 1;
            return r;
        }
        case movearg_First:
        {
            r = num.cbox->move(b, movearg_First);
            assert(r == moveret_OK);
            cursor = 0;
            return r;
        }
        case movearg_Home:
        case movearg_End:
        case movearg_Tab:
        case movearg_Switch:
        case movearg_CtrlSpace:
        {
            if (cursor == 0)
            {
                r = num.cbox->move(b, m);
                assert(r == moveret_OK || r == moveret_End);
                return r;
            }
            else
            {
                assert(cursor == 1);
                r = den.cbox->move(b, m);
                assert(r == moveret_OK || r == moveret_End);
                return r;
            }
        }
        case movearg_Up:
        {
            if (cursor == 0)
            {
                r = num.cbox->move(b, m);
                assert(r == moveret_OK || r == moveret_End);
                return r;
            }
            else
            {
                assert(cursor == 1);
                r = den.cbox->move(b, m);
                assert(r == moveret_OK || r == moveret_End);
                if (r == moveret_OK)
                    return r;
                cursor = 0;
                num.cbox->move(b, movearg_First);
                return moveret_OK;
            }
            
        }
        case movearg_Down:
        {
            if (cursor == 0)
            {
                r = num.cbox->move(b, m);
                assert(r == moveret_OK || r == moveret_End);
                if (r == moveret_OK)
                    return r;
                cursor = 1;
                den.cbox->move(b, movearg_First);
                return moveret_OK;
            }
            else
            {
                assert(cursor == 1);
                r = den.cbox->move(b, m);
                assert(r == moveret_OK || r == moveret_End);
                return r;
            }            
        }
        default:
        {
            assert(false);
            return moveret_OK;
        }
    }
}

insertRet fractionbox::insert(boxbase*&b, insertArg m)
{
    if (cursor == 0)
    {
        return num.cbox->insert(b, m);
    }
    else
    {
        assert(cursor == 1);
        return den.cbox->insert(b, m);
    }
}

removeRet fractionbox::remove(boxbase*&b, removeArg m)
{
    assert(b == nullptr);
    removeRet r;

    if (cursor == 0)
    {
        if (num.cbox->is_selected_placeholder())
        {
            b = nullptr;
            return removeret_Replace;
        }
        else
        {
            r = num.cbox->remove(b, m);
            assert(r == removeret_End || r == removeret_OK);
            if (r == removeret_End)
            {
                b = steal_rowbox(num.cbox, 0,0);
                return removeret_Replace;
            }
            else
            {
                made_into_placeholder(num.cbox);
                return removeret_OK;
            }
        }
    }
    else
    {
        assert(cursor == 1);
        if (den.cbox->is_selected_placeholder())
        {
            b = steal_rowbox(num.cbox, 0,0);
            return removeret_Replace;
        }
        else
        {
            r = den.cbox->remove(b, m);
            assert(r == removeret_End || r == removeret_OK);
            if (r == removeret_End)
            {
                b = steal_rowbox(num.cbox, 0,0);
                return removeret_Replace;
            }
            else
            {
                made_into_placeholder(den.cbox);
                return removeret_OK;
            }
        }
    }
}

ex fractionbox::get_ex()
{
    uex t1(num.cbox->get_ex());
    ex t2 = den.cbox->get_ex();
    return emake_node(gs.sym_sFractionBox.copy(), t1.release(), t2);
}


void fractionbox::get_cursor(aftransform * T)
{
    if (cursor == 0)
    {
        num.cbox->get_cursor(T);
        T->orig_x += num.offx;
        T->orig_y += num.offy;
    }
    else
    {
        assert(cursor == 1);
        den.cbox->get_cursor(T);
        T->orig_x += den.offx;
        T->orig_y += den.offy;        
    }
}

void fractionbox::measure(boxmeasurearg ma)
{
    ma.level++;
    ma.style->start_push_options();
    ma.style->push_opt(opt_FontSize, emake_double(15.0/16*ma.style->get_fontsize()));
    ma.style->finish_push_options();
    num.cbox->measure(ma);
    den.cbox->measure(ma);
    ma.style->pop_options();

    fs = ma.style->get_font();
    fcolor = ma.style->get_fontcolor();

    double frac_over = fontsize_frac_over(fs);
    double frac_under = fontsize_frac_under(fs);
    double frac_center = fontsize_frac_center(fs);

    double extrax = 0.5*frac_center;

    double childmax = std::max(num.cbox->sizex, den.cbox->sizex);
    sizex = 2*extrax + childmax + 2*extrax;

    sizey = num.cbox->sizey + frac_over + frac_under + den.cbox->sizey;
    centery = num.cbox->sizey + frac_center;
    num.offx = 2*extrax + 0.5*(childmax - num.cbox->sizex);
    den.offx = 2*extrax + 0.5*(childmax - den.cbox->sizex);
    num.offy = 0;
    den.offy = num.cbox->sizey + frac_over + frac_under;
}


void align_hline(double & orig_x, double & orig_y, const aftransform* T, double y1, double y2)
{
    if (T->sin_theta == 0)
    {        
        double y1f = (T->orig_y + T->cos_theta*(y1));
        double y2f = (T->orig_y + T->cos_theta*(y2));
        double y1fs = (slong)(y1f + 0.5) - y1f;
        double y2fs = (slong)(y2f + 0.5) - y2f;
        orig_x = T->orig_x;
        orig_y = T->orig_y + 0.75*(std::abs(y1fs) < std::abs(y2fs) ? y1fs : y2fs);
    }
    else if (T->cos_theta == 0)
    {
        double x1f = (T->orig_x - T->sin_theta*y1);
        double x2f = (T->orig_x - T->sin_theta*y2);
        double x1fs = (slong)(x1f + 0.5) - x1f;
        double x2fs = (slong)(x2f + 0.5) - x2f;
        orig_x = T->orig_x + 0.75*(std::abs(x1fs) < std::abs(x2fs) ? x1fs : x2fs);
        orig_y = T->orig_y;
    }
    else
    {
        orig_x = T->orig_x;
        orig_y = T->orig_y;
    }
}

void fractionbox::draw_main(boxdrawarg da)
{
//std::cout << "fractionbox::draw_main " << da.tostring() << std::endl;

    num.cbox->draw_main(boxdrawarg(da, num.offx, num.offy, cursor == 0 ? 0 : DFLAG_IGNORESEL));
    den.cbox->draw_main(boxdrawarg(da, den.offx, den.offy, cursor == 1 ? 0 : DFLAG_IGNORESEL));

    uint32_t color = (da.dflags & DFLAG_SELECTION) ? da.nb->cSelectionForeground :
                     (da.dflags & DFLAG_SCOLOR)    ? da.nb->pallet1[lextype_opinf] :
                                                     fcolor;

    double frac_center = fontsize_frac_center(fs);
    double extrax = 0.5*frac_center;

    double f = fontsize_size(fs);
    double lw = fontsize_line_width(fs);
    double y = da.globy + centery;
    double x = da.globx;
    double e = 0.5*lw;

    double q[5][2];
    q[0][0] = x+extrax;        q[0][1] = y-e;
    q[1][0] = x+extrax;        q[1][1] = y+e;
    q[2][0] = x+sizex-extrax;  q[2][1] = y+e;
    q[3][0] = x+sizex-extrax;  q[3][1] = y-e;
    q[4][0] = x+extrax;        q[4][1] = y-e;


    double orig_x, orig_y;
    align_hline(orig_x, orig_y, da.T, y-e, y+e);

    for (int j = 0; j < 5; j++)
    {
        double u = orig_x + da.T->cos_theta*(q[j][0]) - da.T->sin_theta*(q[j][1]);
        double v = orig_y + da.T->sin_theta*(q[j][0]) + da.T->cos_theta*(q[j][1]);
        q[j][0] = u;
        q[j][1] = v;
    }

    renderPath(&glb_image, q[0], 5, color);
}

void fractionbox::draw_pre(boxdrawarg da)
{
    switch (cursor)
    {
        case 0:
            num.cbox->draw_pre(boxdrawarg(da, num.offx, num.offy));
            return;
        case 1:
            den.cbox->draw_pre(boxdrawarg(da, den.offx, den.offy));
            return;
        default:
            assert(false);
    }
}

void fractionbox::draw_post(boxdrawarg da)
{
    switch (cursor)
    {
        case 0:
            num.cbox->draw_post(boxdrawarg(da, num.offx, num.offy));
            return;
        case 1:
            den.cbox->draw_post(boxdrawarg(da, den.offx, den.offy));
            return;
        default:
            assert(false);
    }
}
