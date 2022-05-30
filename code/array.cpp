#include "globalstate.h"
#include "code.h"
#include "hash.h"
#include "arithmetic.h"
#include "ex_cont.h"

ex dcode_sConstantArray(er e)
{
//std::cout << "dcode_sConstantArray: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sConstantArray.get()));

    if (elength(e) != 2)
        return _handle_message_argx(e, 2);

    er e2 = echild(e,2);
    parray_dims dims;
    if (eis_intnm(e2))
    {
        dims.set_rank(1);
        dims.set_index(0, eintm_get(e2));
    }
    else if (ehas_head_sym(e2, gs.sym_sList.get()))
    {
        ulong r = elength(e2);
        if (r == 0)
            return _handle_message_ilsmn(e, 2);
        dims.set_rank(r);
        for (ulong i = 0; i < r; i++)
        {
            er e2i = echild(e2,i+1);
            if (eis_intnm(e2i))
                dims.set_index(i, eintm_get(e2i));
            else
                return _handle_message_ilsmn(e, 2);
        }
    }
    else if (eis_parray(e2))
    {
        if (!eis_parray_fmpz(e2) || eparray_rank(e2) != 1)
            return _handle_message_ilsmn(e, 2);

        ulong r = eto_parray_fmpz(e2).dimensions.get_index(0);
        assert(r > 0);
        dims.set_rank(r);
        fmpz * d = eto_parray_fmpz(e2).array->data;
        for (ulong i = 0; i < r; i++)
        {
            if (!COEFF_IS_MPZ(d[i]) && d[i] > 0)
                dims.set_index(i, d[i]);
            else
                return _handle_message_ilsmn(e, 2);
        }       
    }
    else
    {
        return _handle_message_ilsmn(e, 2);
    }

    er e1 = echild(e,1);

    ulong size;
    if (dims.rank <= 3 && dims.size_fits_ui() && (size = dims.size()) > 0)
    {
        if (eis_int(e1))
        {
            if (!COEFF_IS_MPZ(*eint_data(e1)))
            {
                uex r(emake_parray_fmpz());
                eto_parray_fmpz(r.get()).dimensions.set(dims);
                fmpz * d = eto_parray_fmpz(r.get()).fit_edit(size);
                for (ulong i = 0; i < size; i++)
                    d[i] = *eint_data(e1);
                return r.release();
            }
        }
        else if (eis_rat(e1))
        {
            if (!COEFF_IS_MPZ(*fmpq_numref(erat_data(e1))) &&
                !COEFF_IS_MPZ(*fmpq_denref(erat_data(e1))))
            {
                uex r(emake_parray_fmpq());
                eto_parray_fmpq(r.get()).dimensions.set(dims);
                fmpq * d = eto_parray_fmpq(r.get()).fit_edit(size);
                for (ulong i = 0; i < size; i++)
                    d[i] = *erat_data(e1);
                return r.release();
            }
        }
        else if (eis_double(e1))
        {
            uex r(emake_parray_double());
            eto_parray_double(r.get()).dimensions.set(dims);
            double * d = eto_parray_double(r.get()).fit_edit(size);
            for (ulong i = 0; i < size; i++)
                d[i] = edouble_number(e1);
            return r.release();
        }

        else if (eis_cmplx_double(e1))
        {
            std::complex<double> z = ecmplx_double_number(e1);
            uex r(emake_parray_cmplx());
            eto_parray_cmplx(r.get()).dimensions.set(dims);
            std::complex<double> * d = eto_parray_cmplx(r.get()).fit_edit(size);
            for (ulong i = 0; i < size; i++)
                d[i] = z;
            return r.release();
        }
    }

    uex r(ecopy(e1));
    for (ulong i = dims.rank; i > 0; i--)
    {
        ulong m = dims.get_index(i - 1);
        uex newr; newr.init_push_backr(gs.sym_sList.get(), m);
        for (ulong j = 0; j < m; j++)
            newr.push_back(r.copy());
        r.setnz(newr.release());
    }
    return r.release();
}


static ex _do_append(er a, er b, er e)
{
    if (eis_node(a))
    {
        size_t n = elength(a);
        uex r; r.init_push_backr(echild(a, 0), n + 1);
        for (size_t i = 1; i <= n; i++)
            r.push_back(ecopychild(a,i));
        r.push_back(ecopy(b));
        return r.release();
    }
    else if (eis_parray(a))
    {
        if (eis_parray_fmpz(a))
        {
            if (eis_int(b))
            {
                uex r(emake_parray_fmpz());
                parray_set(eto_parray_fmpz(r.get()), eto_parray_fmpz(a));
                if (parray_append(eto_parray_fmpz(r.get()), eint_data(b)))
                    return r.release();
            }
            if (eis_parray_fmpz(b))
            {
                uex r(emake_parray_fmpz());
                parray_set(eto_parray_fmpz(r.get()), eto_parray_fmpz(a));
                if (parray_append(eto_parray_fmpz(r.get()), eto_parray_fmpz(b)))
                    return r.release();
            }
            else if (eis_parray_fmpq(b))
            {
                uex r(emake_parray_fmpq());
                parray_set(eto_parray_fmpq(r.get()), eto_parray_fmpz(a));
                if (parray_append(eto_parray_fmpq(r.get()), eto_parray_fmpq(b)))
                    return r.release();
            }
        }
        else if (eis_parray_fmpq(a))
        {
            if (eis_rat(b))
            {
                uex r(emake_parray_fmpq());
                parray_set(eto_parray_fmpq(r.get()), eto_parray_fmpq(a));
                if (parray_append(eto_parray_fmpq(r.get()), erat_data(b)))
                    return r.release();
            }
            else if (eis_parray_fmpz(b))
            {
                uex r(emake_parray_fmpq());
                parray_set(eto_parray_fmpq(r.get()), eto_parray_fmpq(a));
                if (parray_append(eto_parray_fmpq(r.get()), eto_parray_fmpz(b)))
                    return r.release();
            }
            else if (eis_parray_fmpq(b))
            {
                uex r(emake_parray_fmpq());
                parray_set(eto_parray_fmpq(r.get()), eto_parray_fmpq(a));
                if (parray_append(eto_parray_fmpq(r.get()), eto_parray_fmpq(b)))
                    return r.release();
            }
        }
        else if (eis_parray_double(a))
        {
            if (eis_double(b))
            {
                uex r(emake_parray_double());
                parray_set(eto_parray_double(r.get()), eto_parray_double(a));
                if (parray_append(eto_parray_double(r.get()), edouble_number(b)))
                    return r.release();
            }
            else if (eis_parray_double(b))
            {
                uex r(emake_parray_double());
                parray_set(eto_parray_double(r.get()), eto_parray_double(a));
                if (parray_append(eto_parray_double(r.get()), eto_parray_double(b)))
                    return r.release();
            }
            else if (eis_parray_cmplx(b))
            {
                uex r(emake_parray_cmplx());
                parray_set(eto_parray_cmplx(r.get()), eto_parray_double(a));
                if (parray_append(eto_parray_cmplx(r.get()), eto_parray_cmplx(b)))
                    return r.release();
            }
        }
        else if (eis_parray_cmplx(a))
        {
            if (eis_cmplx_double(b))
            {
                uex r(emake_parray_cmplx());
                parray_set(eto_parray_cmplx(r.get()), eto_parray_cmplx(a));
                if (parray_append(eto_parray_cmplx(r.get()), ecmplx_double_number(b)))
                    return r.release();
            }
            else if (eis_parray_cmplx(a) && eis_parray_double(b))
            {
                uex r(emake_parray_cmplx());
                parray_set(eto_parray_cmplx(r.get()), eto_parray_double(a));
                if (parray_append(eto_parray_cmplx(r.get()), eto_parray_double(b)))
                    return r.release();
            }
            else if (eis_parray_cmplx(a) && eis_parray_cmplx(b))
            {
                uex r(emake_parray_cmplx());
                parray_set(eto_parray_cmplx(r.get()), eto_parray_double(a));
                if (parray_append(eto_parray_cmplx(r.get()), eto_parray_cmplx(b)))
                    return r.release();
            }
        }
        er an = eparray_get_normal(a);
        ulong n = elength(an);
        uex r; r.init_push_backr(echild(a, 0), n + 1);
        for (ulong i = 0; i < n; i++)
            r.push_back(ecopychild(an,i+1));
        r.push_back(ecopy(b));
        return r.release();        
    }
    else if (eis_hmap(a))
    {
        uex r(emake_hmap());
        eto_hmap(r.get()).map_data = eto_hmap(a).map_data;
        if (ehas_head_sym(b, gs.sym_sList.get()))
        {
            for (ulong i = 0; i < elength(b); i++)
            {
                er bi = echild(e,i+1);
                if (ehas_head_sym_length(bi, gs.sym_sRule.get(), 2) ||
                    ehas_head_sym_length(bi, gs.sym_sRuleDelayed.get(), 2))
                {
                    ehmap_assign(r.get(), ecopychild(bi,1), ecopychild(bi,2));
                }
                else
                {
                    _gen_message(esymbolic_head(e), "invdt", "The argument `1` is not a rule or a list of rules.", emake_cint(1), ecopy(b));
                    return ecopy(e);
                }
            }
        }
        else if (ehas_head_sym_length(b, gs.sym_sRule.get(), 2) ||
                 ehas_head_sym_length(b, gs.sym_sRuleDelayed.get(), 2))
        {
            ehmap_assign(r.get(), ecopychild(b,1), ecopychild(b,2));
        }
        else
        {
            _gen_message(esymbolic_head(e), "invdt", "The argument `1` is not a rule or a list of rules.", emake_cint(1), ecopy(b));
            return ecopy(e);
        }
        return r.release();
    }
    else
    {
        _gen_message(esymbolic_head(e), "normal", NULL, emake_cint(1), ecopy(e));
        return ecopy(e);
    }
}


ex scode_sAppend(er e)
{
//std::cout << "scode_sAssociation: " << ex_tostring_full(e) << std::endl;
    assert(eis_node(e));

    er h = echild(e,0);
    if (!ehas_head_sym_length(h, gs.sym_sList.get(), 1))
        return ecopy(e);

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    return _do_append(echild(h,1), echild(e,1), e);
}

ex dcode_sAppend(er e)
{
    assert(ehas_head_sym(e, gs.sym_sAppend.get()));

    if (elength(e) == 2)
    {
        return _do_append(echild(e,1), echild(e,2), e);
    }
    else if (elength(e) == 1)
    {
        return ecopy(e);
    }
    else
    {
        return _handle_message_argt(e, (1 << 0) + (2 << 8));
    }
}

// return < 0 for error, other the number of steps
slong _parse_range(
    uex& lower, uex& upper, uex& step,
    er it, slong ioff, slong itn)
{
    if (itn == 1)
    {
        lower.setz(emake_cint(1));
        upper.setz(ecopychild(it, ioff + 1));
        step.setz(emake_cint(1));
    }
    else if (itn == 2)
    {
        lower.setz(ecopychild(it, ioff + 1));
        upper.setz(ecopychild(it, ioff + 2));
        step.setz(emake_cint(1));

    }
    else if (itn == 3)
    {
        lower.setz(ecopychild(it, ioff + 1));
        upper.setz(ecopychild(it, ioff + 2));
        step.setz(ecopychild(it, ioff + 3));
    }
    else
    {
        return -1;
    }

    uex scount(ex_subx(upper.copy(), lower.copy()));
    scount.setz(ex_divx(scount.release(), step.copy()));
    scount.setz(ex_floor(scount.release()));
    if (eis_int(scount.get()))
    {
        if (fmpz_sgn(scount.int_data()) < 0)
            return 0;

        if (eis_intsm(scount.get()))
            return 1 + eintsm_get(scount.get());
    }
    return -1;
}


ex dcode_sRange(er e)
{
//std::cout << "dcode_sRange: " << ex_tostring_full(e) << std::endl;
    assert(eis_node(e));

    if (elength(e) < 1 || elength(e) > 3)
        return _handle_message_argb(e, 1 + (3 << 8));

    uex lower, upper, step;
    slong n = _parse_range(lower, upper, step, e, 0, elength(e));
    if (n < 0)
    {
        _gen_message(echild(e,0), "range", "Range specification `1` does not have appropriate bounds.", ecopy(e));
        return ecopy(e);
    }
    else if (n == 0)
    {
        return gs.const_empty_list.copy();
    }

    uex r; r.init_push_backr(gs.sym_sList.get(), n);
    for (slong i = 0; i < n; i++)
        r.push_back(ex_addx(lower.copy(), ex_mul_six(step.copy(), i)));
    r.setnz(econvert_parray_shallow(r.get()));
    return r.release();
}
