#include <stack>
#include "boxes.h"
#include "ex_parse.h"
#include "notebook.h"
#include "box_convert.h"

void sqrtbox::print(size_t indent, double offx, double offy)
{
    for (size_t i = 0; i < indent; i++)
        printf("    ");
    printf("(%f,%f) sqrt: size (%f,%f:%f)\n", offx, offy, sizex, sizey, centery);

    inside.cbox->print(indent + 1, inside.offx, inside.offy);
}

visitRet sqrtbox::visit(visitArg m)
{
    inside.cbox->visit(m);
    return visitret_OK;
}

boxbase * sqrtbox::copy()
{
    rowbox * newinside = dynamic_cast<rowbox*>(inside.cbox->copy());
    sqrtbox * r = new sqrtbox(newinside);
    r->inside.offx = inside.offx;
    r->inside.offy = inside.offy;
    return r;    
}

void sqrtbox::key_copy(boxbase*&b)
{
    inside.cbox->key_copy(b);
}

void sqrtbox::key_paste(boxbase*&b)
{
    inside.cbox->key_paste(b);
}

void sqrtbox::insert_char(int32_t c)
{
    inside.cbox->insert_char(c);
    return;
}

moveRet sqrtbox::move(boxbase*&b, moveArg m)
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

insertRet sqrtbox::insert(boxbase*&b, insertArg m)
{
    return inside.cbox->insert(b, m);
}

removeRet sqrtbox::remove(boxbase*&b, removeArg m)
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

ex sqrtbox::get_ex()
{
    ex t = inside.cbox->get_ex();
    return emake_node(gs.sym_sSqrtBox.copy(), t);
}


void sqrtbox::get_cursor(aftransform * T)
{
    inside.cbox->get_cursor(T);
    T->orig_x += inside.offx;
    T->orig_y += inside.offy;    
    return;
}


void sqrtbox::measure(boxmeasurearg ma)
{
    ma.level++;
    inside.cbox->measure(ma);

    fs = ma.style->get_font();
    fcolor = ma.style->get_fontcolor();

    double padup    = fontsize_sqrt_padup(fs);
    double f = fontsize_size(fs);
    double paddown  = f*0.125;
    double padleft  = f*1.625;
    double padright = f*0.375;

    inside.offx = padleft;
    inside.offy = padup;
    sizex = padleft + inside.cbox->sizex + padright;
    sizey = padup + inside.cbox->sizey + paddown;
    centery = padup + inside.cbox->centery;
}

/*
lineleft[{x1_,y1_},{x2_,y2_},\[Epsilon]_]:=
    {y1-y2,x2-x1,x1 y2-x2 y1+\[Epsilon] Sqrt[(x1-x2)^2+(y1-y2)^2]};
*/
void lineleft(
    double & a, double & b, double & c,
    double * p1,
    double * p2,
    double epsilon)
{
    double x1 = p1[0];
    double y1 = p1[1];
    double x2 = p2[0];
    double y2 = p2[1];
    double A = y1 - y2;
    double B = x2 - x1;
    a = A;
    b = B;
    c = x1*y2 - x2*y1 - epsilon*sqrt(A*A + B*B);
}

/*
    Solve[{a1 x + b1 y + c1 == 0, a2 x + b2 y + c2 == 0}, {x, y}]
*/
void intersection(
    double * p,
    double a1, double b1, double c1,
    double a2, double b2, double c2)
{
    p[0] = (b2*c1 - b1*c2)/(a2*b1 - a1*b2);
    p[1] = (a1*c2 - a2*c1)/(a2*b1 - a1*b2);
}

/*
    Solve[{
        lineleft[{x1,y1},{x2,y2},\[Epsilon]].{x,y,1}==0,
        (y-y2)/(x-x2)==-1/((y1-y2)/(x1-x2))
    },{x,y}]//Simplify//CForm
*/
void turnleft(
    double * p,
    double * p1,
    double * p2,
    double epsilon)
{
    double x1 = p1[0];
    double y1 = p1[1];
    double x2 = p2[0];
    double y2 = p2[1];
    double d = (x2-x1)*(x2-x1) + (y1-y2)*(y1-y2);
    double s = sqrt(d);

    p[0] = x2 + epsilon/s*(y1 - y2);
    p[1] = y2 + epsilon/s*(x2 - x1);
}


void sqrtbox::draw_main(boxdrawarg da)
{
    inside.cbox->draw_main(boxdrawarg(da, inside.offx, inside.offy));

    uint32_t color = (da.dflags & DFLAG_SELECTION) ? da.nb->cSelectionForeground :
                     (da.dflags & DFLAG_SCOLOR)    ? da.nb->pallet1[lextype_opinf] :
                                                     fcolor;

    double x = inside.cbox->sizex;
    double y = inside.cbox->sizey;

    double e1 = 0.5*fontsize_line_width(fs);
    double e2 = 0.5*fontsize_sqrt_width(fs);

    double left = inside.offx;
    double up = inside.offy;
    double down = sizey - inside.offy - y;
    double right = sizex - left - x;


    double p[5][2] = {
        {-1.5000/1.625*left   , 0.875*y - 0.875/1.625*left},
        {-1.1875/1.625*left   , 0.875*y - 1.000/1.625*left},
        {-0.7500/1.625*left   , y},
        {-0.3750/1.625*left   , -0.75*up},
        {x + 0.5*right        , -0.75*up},
    };
    for (int i = 0; i <= 4; i++)
    {
        p[i][0] += da.globx + left;
        p[i][1] += da.globy + up;
    }

    double q[15][2];
    double a1,b1,c1, a2,b2,c2;

    // left radical
    turnleft(q[0], p[1], p[0], -e1);

    turnleft(q[1], p[1], p[0], +e1);

    lineleft(a1,b1,c1, p[0], p[1], -e1);
    lineleft(a2,b2,c2, p[1], p[2], -e2);
    intersection(q[2], a1,b1,c1, a2,b2,c2);

    lineleft(a1,b1,c1, p[1], p[2], -e2);
    lineleft(a2,b2,c2, p[2], p[3], -e1);
    intersection(q[3], a1,b1,c1, a2,b2,c2);

    lineleft(a1,b1,c1, p[2], p[3], -e1);
    lineleft(a2,b2,c2, p[3], p[4], -e1);
    intersection(q[4], a1,b1,c1, a2,b2,c2);

    lineleft(a1,b1,c1, p[3], p[4], +e1);
    lineleft(a2,b2,c2, p[3], p[2], -e1);
    intersection(q[5], a1,b1,c1, a2,b2,c2);

    lineleft(a1,b1,c1, p[3], p[2], -e1);
    intersection(q[6], a1,b1,c1, 0,-1,p[2][1]+(e1+e2)/2);

    lineleft(a1,b1,c1, p[2], p[1], -e2);
    intersection(q[7], a1,b1,c1, 0,-1,p[2][1]+(e1+e2)/2);

    lineleft(a1,b1,c1, p[2], p[1], -e2);
    lineleft(a2,b2,c2, p[1], p[0], -e1);
    intersection(q[8], a1,b1,c1, a2,b2,c2);

    q[9][0]=q[0][0]; q[9][1]=q[0][1];

    // top bar
    turnleft(q[10], p[3], p[4], -e1);

    turnleft(q[11], p[3], p[4], +e1);

    lineleft(a1,b1,c1, p[4], p[3], -e1);
    lineleft(a2,b2,c2, p[3], p[2], 0);
    intersection(q[12], a1,b1,c1, a2,b2,c2);

    q[13][0] = q[4][0]; q[13][1] = q[4][1];

    q[14][0] = q[10][0]; q[14][1] = q[10][1];

    double orig_x, orig_y;
    align_hline(orig_x, orig_y, da.T, da.globy + 0.25*up - e1, da.globy + 0.25*up + e1);

    for (int j = 0; j < 15; j++)
    {
        double u = orig_x + da.T->cos_theta*(q[j][0]) - da.T->sin_theta*(q[j][1]);
        double v = orig_y + da.T->sin_theta*(q[j][0]) + da.T->cos_theta*(q[j][1]);
        q[j][0] = u;
        q[j][1] = v;
    }

    renderPath(&glb_image, q[0], 10, color);
    renderPath(&glb_image, q[10], 5, color);
}

void sqrtbox::draw_pre(boxdrawarg da)
{
    inside.cbox->draw_pre(boxdrawarg(da, inside.offx, inside.offy));
}

void sqrtbox::draw_post(boxdrawarg da)
{
    inside.cbox->draw_post(boxdrawarg(da, inside.offx, inside.offy));
}
