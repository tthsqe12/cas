#include <stack>
#include "boxes.h"
#include "ex_parse.h"
#include "notebook.h"
#include "box_convert.h"

void tagbox::print(size_t indent, double offx, double offy)
{
    for (size_t i = 0; i < indent; i++)
        printf("    ");
    printf("(%f,%f) tag: size (%f,%f:%f)\n", offx, offy, sizex, sizey, centery);

    inside.cbox->print(indent + 1, inside.offx, inside.offy);
}

visitRet tagbox::visit(visitArg m)
{
    inside.cbox->visit(m);
    return visitret_OK;
}

boxbase * tagbox::copy()
{
    rowbox * newinside = dynamic_cast<rowbox*>(inside.cbox->copy());
    tagbox * r = new tagbox(newinside, tag.copy());
    r->inside.offx = inside.offx;
    r->inside.offy = inside.offy;
    return r;    
}

void tagbox::key_copy(boxbase*&b)
{
    inside.cbox->key_copy(b);
}

void tagbox::key_paste(boxbase*&b)
{
    inside.cbox->key_paste(b);
}

void tagbox::insert_char(int32_t c)
{
    inside.cbox->insert_char(c);
    return;
}

moveRet tagbox::move(boxbase*&b, moveArg m)
{
    moveRet r;
    switch (m)
    {
        case movearg_Left:
        case movearg_ShiftLeft:
        case movearg_Right:
        case movearg_ShiftRight:
        case movearg_Home:
        case movearg_End:
        case movearg_Up:
        case movearg_Down:
        {
            r = inside.cbox->move(b, m);
            assert(r == moveret_OK || r == moveret_End);
            return r;
        }
        case movearg_Last:
        case movearg_First:
        {
            r = inside.cbox->move(b, m);
            assert(r == moveret_OK);
            return r;
        }
        case movearg_CtrlSpace:
        {
            r = inside.cbox->move(b, m);
            assert(r == moveret_OK || r == moveret_End);
            return r;            
        }
        case movearg_Switch:
        {
            inside.cbox->move(b, m);
            return moveret_OK;
        }
        default:
        {
            assert(false);
            return moveret_OK;
        }
    }
}

insertRet tagbox::insert(boxbase*&b, insertArg m)
{
    return inside.cbox->insert(b, m);
}

removeRet tagbox::remove(boxbase*&b, removeArg m)
{
    assert(b == nullptr);
    removeRet r;

    if (inside.cbox->is_selected_placeholder())
    {
        b = nullptr;
        return removeret_Replace;
    }
    else
    {
        r = inside.cbox->remove(b, m);
        assert(r == removeret_End || r == removeret_OK);
        if (r == removeret_End)
        {
            b = steal_rowbox(inside.cbox, 0,0);
            return removeret_Replace;
        }
        else
        {
            made_into_placeholder(inside.cbox);
            return removeret_OK;
        }
    }
}

ex tagbox::get_ex()
{
    ex t = inside.cbox->get_ex();
    return emake_node(gs.sym_sTagBox.copy(), t, tag.copy());
}


void tagbox::get_cursor(aftransform * T)
{
    inside.cbox->get_cursor(T);
    T->orig_x += inside.offx;
    T->orig_y += inside.offy;    
    return;
}


void tagbox::measure(boxmeasurearg ma)
{
    ma.level++;
    inside.cbox->measure(ma);

    double padup    = 0;
    double paddown  = 0;
    double padleft  = 0;
    double padright = 0;

    inside.offx = padleft;
    inside.offy = padup;
    sizex = padleft + inside.cbox->sizex + padright;
    sizey = padup + inside.cbox->sizey + paddown;
    centery = padup + inside.cbox->centery;
}


void tagbox::draw_main(boxdrawarg da)
{
    inside.cbox->draw_main(boxdrawarg(da, inside.offx, inside.offy));
}

void tagbox::draw_pre(boxdrawarg da)
{
    inside.cbox->draw_pre(boxdrawarg(da, inside.offx, inside.offy));
}

void tagbox::draw_post(boxdrawarg da)
{
    inside.cbox->draw_post(boxdrawarg(da, inside.offx, inside.offy));
}
