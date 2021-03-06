#include "globalstate.h"
#include "code.h"
#include "ex_print.h"
#include "eval.h"
#include "ex_cont.h"

static ex det(er m, std::vector<size_t> & row, std::vector<size_t> & col)
{
    assert(row.size() >= 1);
    assert(row.size() == col.size());
    size_t n = row.size();

    if (n == 1)
    {
        return ecopychild(m, row[0], col[0]);
    }

    std::vector<size_t> nrow, ncol;
    for (size_t i = 1; i < n; i++)
    {
        nrow.push_back(row[i]);
        ncol.push_back(col[i]);
    }

    uex r; r.init_push_backr(gs.sym_sPlus.get(), n);
    for (size_t i = 0; i < n; i++)
    {
        ex s = det(m, nrow, ncol);
        ex c = ecopychild(m, row[0], col[i]);
        s = (i%2) == 0 ? ex_mulx(c, s) : ex_mulx(emake_cint(-1), c, s);
        r.push_back(s);
        if (i + 1 < n)
            ncol[i] = col[i];
    }

    return ex_canonicalize_plus(r.release());
}

ex dcode_sInverse(er e)
{
//std::cout << "dcode_sInverse: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sInverse.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    er m = echild(e,1);
    
    if (!ehas_head_sym(m, gs.sym_sList.get()))
        return ecopy(e);

    size_t n = elength(m);
    if (n == 0)
        return ecopy(e);

    for (size_t i = 1; i <= n; i++)
    {
        if (!ehas_head_sym_length(echild(m,i), gs.sym_sList.get(), n))
            return ecopy(e);
    }

    std::vector<size_t> row, col;
    for (size_t i = 1; i <= n; i++)
    {
        row.push_back(i);
        col.push_back(i);
    }
    uex d(det(m, row, col));

//std::cout << "det: " << ex_tostring(d.get()) << std::endl;

    uex im; im.init_push_backr(gs.sym_sList.get(), n);
    for (size_t i = 1; i<= n; i++)
    {
        uex irow; irow.init_push_backr(gs.sym_sList.get(), n);
        for (size_t j = 1; j <= n; j++)
        {
            row.clear();
            col.clear();
            for (size_t k = 1; k <= n; k++)
            {
                if (k != j)
                    row.push_back(k);
                if (k != i)
                    col.push_back(k);
            }

            ex c = det(m, row, col);
            if (((i + j) % 2) != 0)
            {
                c = ex_mulx(emake_cint(-1), c);
            }
            irow.push_back(ex_divx(c, d.copy()));
        }
        im.push_back(irow.release());
    }

    return im.release();
}


static ex _unpack(er e, ulong depth, ulong maxdepth)
{
	if (depth > maxdepth)
		return ecopy(e);

	uex r(ecopy(e));
	if (eis_parray(e))
		r.setnz(ecopy(eparray_get_normal(e)));

	if (eis_node(r.get()))
	{
		ulong n = elength(r.get());
		uex f; f.init_push_backx(_unpack(r.child(0), depth + 1, maxdepth), n);
		for (ulong i = 1; i <= n; i++)
			f.push_back(_unpack(r.child(i), depth + 1, maxdepth));
		r.setnz(f.release());
	}

	return r.release();
}


static ex dot(er aa, er bb)
{
    uex A(ecopy(aa)), B(ecopy(bb));
    if (eis_parray(aa))
        A.setnz(_unpack(aa, 1, -UWORD(1)));
    if (eis_parray(bb))
        B.setnz(_unpack(bb, 1, -UWORD(1)));

    er a = A.get();
    er b = B.get();

    if (!ehas_head_sym(a, gs.sym_sList.get()) ||
        !ehas_head_sym(b, gs.sym_sList.get()))
    {
        return nullptr;
    }

    size_t listcount = 0;
    for (size_t i = 1; i <= elength(a); i++)
    {
        listcount += ehas_head_sym(echild(a,i), gs.sym_sList.get());
    }

    if (listcount == 0)
    {
        if (elength(b) != elength(a))
            return nullptr;

        uex r; r.init_push_backr(gs.sym_sPlus.get(), elength(a));
        for (size_t i = 1; i <= elength(a); i++)
            r.push_back(eval(emake_node(gs.sym_sTimes.copy(), ecopychild(a,i), ecopychild(b,i))));

        return eval(r.release());
    }
    else if (listcount == elength(a))
    {
        uex r; r.init_push_backr(gs.sym_sList.get(), elength(a));
        for (size_t i = 1; i <= elength(a); i++)
        {
            ex s = dot(echild(a,i), b);
            if (s == nullptr)
                return nullptr;

            r.push_back(s);
        }
        return r.release();
    }
    else
    {
        return nullptr;
    }
}

ex dcode_sDot(er e)
{
//std::cout << "dcode_sDot: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sDot.get()));

    if (elength(e) == 0)
    {
        return ecopy(e);
    }

    uex r(ecopychild(e,1));
    for (ulong i = 1; i < elength(e); i++)
    {
        r.setnz(dot(r.get(), echild(e,i+1)));
        if (r.get() == nullptr)
        {
            return ecopy(e);
        }
    }

    return r.release();
}

ex dcode_sIdentityMatrix(er e)
{
//std::cout << "dcode_sIdentityMatrix: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sIdentityMatrix.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    if (!eis_intpm(echild(e,1)))
    {
		_gen_message(echild(e,0), "intpm", NULL, ecopy(e), emake_cint(1));
        return ecopy(e);
    }

    slong n = eintm_get(echild(e,1));
	slong msize;
	if (mul_si_checked(msize, n, n))
	{
        _gen_message(echild(e,0), "toomany", NULL, ecopy(e));
        return ecopy(e);
	}

    uex m(emake_parray_fmpz());
    eto_parray_fmpz(m.get()).dimensions.set_rank2(n, n);
	fmpz * M = eto_parray_fmpz(m.get()).fit_edit(msize);
    for (slong i = 0; i < n; i++)
		fmpz_one(M + i*n + i);
    return m.release();
}
