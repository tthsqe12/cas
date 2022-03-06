#include <stack>
#include "boxes.h"
#include "ex_parse.h"
#include "notebook.h"
#include "box_convert.h"

void subscriptbox::print(size_t indent, double offx, double offy)
{
    for (size_t i = 0; i < indent; i++)
        printf("    ");

    printf("(%f,%f) subscript: size(%f,%f:%f)\n",
             offx, offy, sizex, sizey, sizey);

    sub.cbox->print(indent + 1, sub.offx, sub.offy);
}

visitRet subscriptbox::visit(visitArg m)
{
    sub.cbox->visit(m);
    return visitret_OK;
}

boxbase * subscriptbox::copy()
{
    rowbox * newsub = dynamic_cast<rowbox*>(sub.cbox->copy());
    subscriptbox * r = new subscriptbox(newsub);
    r->sub.offx = sub.offx;
    r->sub.offy = sub.offy;
    return r;    
}

void subscriptbox::key_copy(boxbase*&b)
{
    sub.cbox->key_copy(b);
}

void subscriptbox::key_paste(boxbase*&b)
{
    sub.cbox->key_paste(b);
}

void subscriptbox::insert_char(int32_t c)
{
    sub.cbox->insert_char(c);
}

moveRet subscriptbox::move(boxbase*&b, moveArg m)
{
    moveRet r;
    switch (m)
    {
        case movearg_Left:
        case movearg_ShiftLeft:
        case movearg_Right:
        case movearg_ShiftRight:
        {
            r = sub.cbox->move(b, m);
            assert(r == moveret_OK || r == moveret_End);
            return r;
        }
        case movearg_Last:
        case movearg_First:
        {
            r = sub.cbox->move(b, m);
            assert(r == moveret_OK);
            return r;
        }
        case movearg_Home:
        case movearg_End:
        case movearg_Tab:
        case movearg_CtrlSpace:
        {
            r = sub.cbox->move(b, m);
            assert(r == moveret_OK || r == moveret_End);
            return r;
        }
        case movearg_Switch:
        {
            rowbox * newsub = steal_rowbox(sub.cbox, 0,0);
            rowbox* newsuper = new rowbox(2, 0,1);
            newsuper->child[0].cibox = iboximm_make(CHAR_Placeholder);
            newsuper->child[1].cibox.ptr = new nullbox();
            b = new subsuperscriptbox(newsub, newsuper, 1);
            return moveret_Replace;
        }
        case movearg_Up:
        case movearg_Down:
        {
            r = sub.cbox->move(b, m);
            assert(r == moveret_OK || r == moveret_End);
            return r;
        }
        default:
        {
            assert(false);
            return moveret_OK;
        }
    }
}

insertRet subscriptbox::insert(boxbase*&b, insertArg m)
{
    return sub.cbox->insert(b, m);
}

removeRet subscriptbox::remove(boxbase*&b, removeArg m)
{
    removeRet r;
    if (sub.cbox->is_selected_placeholder())
    {
        b = nullptr;
        return removeret_Replace;
    }
    else
    {
        r = sub.cbox->remove(b, m);
        assert(r == removeret_End || r == removeret_OK);
        if (made_into_placeholder(sub.cbox))
        {
            return removeret_OK;
        }
        if (r == removeret_End)
        {
            b = steal_rowbox(sub.cbox, 0,0);
            return removeret_Replace;
        }
        else
        {
            return removeret_OK;
        }
    }
}

ex subscriptbox::get_ex()
{
    ex t = sub.cbox->get_ex();
    return emake_node(gs.sym_sSubscriptBox.copy(), t);
}


void subscriptbox::get_cursor(aftransform * T)
{
    sub.cbox->get_cursor(T);
    T->orig_x += sub.offx;
    T->orig_y += sub.offy;
}

void subscriptbox::measure(boxmeasurearg ma)
{
    ma.style->new_options();
    ma.style->push_opt(opt_FontSize, emake_double(0.75*ma.style->get_fontsize()));
    ma.style->push_opt(opt_AutoSpacing, gs.sym_sFalse.copy());
    ma.style->init_const1();
    ma.level++;
    sub.cbox->measure(ma);
    ma.style->pop_options();

    sub.offx = 0;
    sizex = sub.cbox->sizex;
    centery = sub.cbox->sizey;
    sizey = centery + ma.style->get_default_centery();
    sub.offy = 0;
}

void subscriptbox::draw_pre(boxdrawarg da)
{
    sub.cbox->draw_pre(boxdrawarg(da, sub.offx, sub.offy));
}

void subscriptbox::draw_main(boxdrawarg da)
{
    sub.cbox->draw_main(boxdrawarg(da, sub.offx, sub.offy));
}

void subscriptbox::draw_post(boxdrawarg da)
{
    sub.cbox->draw_post(boxdrawarg(da, sub.offx, sub.offy));
}
