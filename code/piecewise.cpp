#include "globalstate.h"
#include "code.h"

ex dcode_sDivisible(er e)
{
//std::cout << "dcode_sDivisible: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sDivisible.get()));

    if (elength(e) < 2)
        return _handle_message_argm(e, 2);

    return ecopy(e);
}

ex dcode_sSign(er e)
{
//std::cout << "dcode_sSign: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sSign.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    er x = echild(e,1);

    switch (etype(x))
    {
        case ETYPE_INT:
            return emake_cint(fmpz_sgn(eint_data(x)));
        case ETYPE_RAT:
            return emake_cint(fmpz_sgn(fmpq_numref(erat_data(x))));
        case ETYPE_DOUBLE:
            return emake_cint(edouble_number(x) > 0 ? 1 : edouble_number(x) < 0 ? -1 : 0);
        case ETYPE_REAL:
            if (arb_is_positive(ereal_data(x)))
                return emake_cint(1);
            else if (arb_is_negative(ereal_data(x)))
                return emake_cint(-1);
            else if (arb_is_zero(ereal_data(x)))
                return emake_cint(0);
            else
                return gs.const_zero_pm_one.copy();
        default:
            return ecopy(e);
    }
}

ex dcode_sAbs(er e)
{
//std::cout << "dcode_sAbs: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sAbs.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    er x = echild(e,1);

    switch (etype(x))
    {
        case ETYPE_INT:
        {
            if (fmpz_sgn(eint_data(x)) >= 0)
            {
                return ecopy(x);
            }
            else
            {
                ex z = emake_int();
                fmpz_neg(eint_data(z), eint_data(x));
                return efix_int(z);
            }
        }
        case ETYPE_RAT:
        {
            if (fmpq_sgn(erat_data(x)) >= 0)
            {
                return ecopy(x);
            }
            else
            {
                ex z = emake_rat();
                fmpq_neg(erat_data(z), erat_data(x));
                return efix_rat(z);
            }
        }
        case ETYPE_DOUBLE:
        {
            if (edouble_number(x) >= 0)
            {
                return ecopy(x);
            }
            else
            {
                return emake_double(-edouble_number(x));
            }
        }
        case ETYPE_REAL:
        {
            if (!arb_contains_negative(ereal_data(x)))
            {
                return ecopy(x);
            }
            else
            {
                ex z = emake_real();
                arb_abs(ereal_data(z), ereal_data(x));
                return z;
            }
        }
        default:
        {
            return ecopy(e);
        }
    }
}


ex dcode_sFloor(er e)
{
//std::cout << "dcode_sFloor: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sFloor.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    er x = echild(e,1);

    switch (etype(x))
    {
        case ETYPE_INT:
        {
            return ecopy(x);
        }
        case ETYPE_RAT:
        {
            ex z = emake_int();
            fmpz_fdiv_q(eint_data(z), fmpq_numref(erat_data(x)), fmpq_denref(erat_data(x)));
            return efix_int(z);
        }
        case ETYPE_DOUBLE:
        {
            ex z = emake_int();
            nice_ball y;
            y.set_double(edouble_number(x));
            if (y.bad)
                return gs.const_indeterminate.copy();
            fmpz_mul_si(y.center, y.center, y.sign);
            if (y.exp > 0)
                fmpz_fdiv_q_2exp(eint_data(z), y.center, y.exp);
            else
                fmpz_mul_2exp(eint_data(z), y.center, -y.exp);
            return z;
        }
        case ETYPE_REAL:
        {
            xarb_t y;
            ex z = emake_int();
            slong p = ereal_number(x).wprec();
            arb_floor(y.data, ereal_data(x), p + EXTRA_PRECISION_BASIC);
            if (arb_get_unique_fmpz(eint_data(z), y.data))
                return z;
            eclear(z);
            return emake_real_move(y);
        }
        default:
        {
            return ecopy(e);
        }
    }
}

ex dcode_sCeiling(er e)
{
//std::cout << "dcode_sCeiling: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sCeiling.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    er x = echild(e,1);

    switch (etype(x))
    {
        case ETYPE_INT:
        {
            return ecopy(x);
        }
        case ETYPE_RAT:
        {
            ex z = emake_int();
            fmpz_cdiv_q(eint_data(z), fmpq_numref(erat_data(x)), fmpq_denref(erat_data(x)));
            return efix_int(z);
        }
        case ETYPE_DOUBLE:
        {
            ex z = emake_int();
            nice_ball y;
            y.set_double(edouble_number(x));
            if (y.bad)
                return gs.const_indeterminate.copy();
            fmpz_mul_si(y.center, y.center, y.sign);
            if (y.exp > 0)
                fmpz_cdiv_q_2exp(eint_data(z), y.center, y.exp);
            else
                fmpz_mul_2exp(eint_data(z), y.center, -y.exp);
            return z;
        }
        case ETYPE_REAL:
        {
            xarb_t y;
            ex z = emake_int();
            slong p = ereal_number(x).wprec();
            arb_ceil(y.data, ereal_data(x), p + EXTRA_PRECISION_BASIC);
            if (arb_get_unique_fmpz(eint_data(z), y.data))
                return z;
            eclear(z);
            return emake_real_move(y);
        }
        default:
        {
            return ecopy(e);
        }
    }
}


ex dcode_sIntegerPart(er e)
{
//std::cout << "dcode_sIntegerPart: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sIntegerPart.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    er x = echild(e,1);

    switch (etype(x))
    {
        case ETYPE_INT:
        {
            return ecopy(x);
        }
        case ETYPE_RAT:
        {
            ex z = emake_int();
            fmpz_tdiv_q(eint_data(z), fmpq_numref(erat_data(x)), fmpq_denref(erat_data(x)));
            return efix_int(z);
        }
        case ETYPE_DOUBLE:
        {
            ex z = emake_int();
            nice_ball y;
            y.set_double(edouble_number(x));
            if (y.bad)
                return gs.const_indeterminate.copy();
            fmpz_mul_si(y.center, y.center, y.sign);
            if (y.exp > 0)
                fmpz_tdiv_q_2exp(eint_data(z), y.center, y.exp);
            else
                fmpz_mul_2exp(eint_data(z), y.center, -y.exp);
            return z;
        }
        case ETYPE_REAL:
        {
            xarb_t y;
            ex z = emake_int();
            slong p = ereal_number(x).wprec();
            if (arb_is_nonnegative(ereal_data(x)))
                arb_floor(y.data, ereal_data(x), p + EXTRA_PRECISION_BASIC);
            else if (arb_is_nonpositive(ereal_data(x)))
                arb_ceil(y.data, ereal_data(x), p + EXTRA_PRECISION_BASIC);
            else
                arb_set(y.data, ereal_data(x));
            if (arb_get_unique_fmpz(eint_data(z), y.data))
                return z;
            eclear(z);
            return emake_real_move(y);
        }
        default:
        {
            return ecopy(e);
        }
    }
}

ex dcode_sRound(er e)
{
//std::cout << "dcode_sRound: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sRound.get()));

    return ecopy(e);
}

ex dcode_sFractionalPart(er e)
{
//std::cout << "dcode_sFractionalPart: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sFractionalPart.get()));

    return ecopy(e);
}


ex dcode_sMod(er e)
{
//std::cout << "dcode_sMod: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sMod.get()));

    if (elength(e) == 2)
    {
        er x = echild(e,1);
        er y = echild(e,2);
        if (eis_int(x) && eis_int(y))
        {
            if (unlikely(fmpz_is_zero(eint_data(y))))
            {
                _gen_message(echild(e,0), "indet", NULL, ecopy(e));
                return gs.const_indeterminate.copy();
            }
            ex z = emake_int();
            fmpz_fdiv_r(eint_data(z), eint_data(x), eint_data(y));
            return efix_int(z);
        }
        else if (eis_int(x) && eis_rat(y))
        {
            return ecopy(e);
        }
        else if (eis_rat(x) && eis_int(y))
        {
            return ecopy(e);
        }
        else if (eis_rat(x) && eis_rat(y))
        {
            return ecopy(e);
        }
        else
        {
            return ecopy(e);
        }
    }
    else
    {
        return ecopy(e);
    }
    return ecopy(e);
}

ex dcode_sQuotient(er e)
{
//std::cout << "dcode_sQuotient: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sQuotient.get()));

    if (elength(e) == 2)
    {
        er x = echild(e,1);
        er y = echild(e,2);
        if (eis_int(x) && eis_int(y))
        {
            if (unlikely(fmpz_is_zero(eint_data(y))))
            {
                if (fmpz_is_zero(eint_data(x)))
                {
                    _gen_message(echild(e,0), "indet", NULL, ecopy(e));
                    return gs.const_indeterminate.copy();
                }
                else
                {
                    _gen_message(echild(e,0), "infy", NULL, ecopy(e));
                    return gs.const_complexinfinity.copy();
                }
            }
            ex z = emake_int();
            fmpz_fdiv_q(eint_data(z), eint_data(x), eint_data(y));
            return efix_int(z);
        }
        else if (eis_int(x) && eis_rat(y))
        {
            return ecopy(e);
        }
        else if (eis_rat(x) && eis_int(y))
        {
            return ecopy(e);
        }
        else if (eis_rat(x) && eis_rat(y))
        {
            return ecopy(e);
        }
        else
        {
            return ecopy(e);
        }
    }
    else
    {
            return ecopy(e);
    }
    return ecopy(e);
}

ex dcode_sQuotientRemainder(er e)
{
//std::cout << "dcode_sQuotient: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sQuotientRemainder.get()));

    if (elength(e) == 2)
    {
        er x = echild(e,1);
        er y = echild(e,2);
        if (eis_int(x) && eis_int(y))
        {
            if (unlikely(fmpz_is_zero(eint_data(y))))
            {
                _gen_message(echild(e,0), "divz", NULL, ecopy(y), ecopy(e));
                return ecopy(e);
            }
            ex q = emake_int();
            ex r = emake_int();
            fmpz_fdiv_qr(eint_data(q), eint_data(r), eint_data(x), eint_data(y));
            return emake_list(efix_int(q), efix_int(r));
        }
        else if (eis_int(x) && eis_rat(y))
        {
            return ecopy(e);
        }
        else if (eis_rat(x) && eis_int(y))
        {
            return ecopy(e);
        }
        else if (eis_rat(x) && eis_rat(y))
        {
            return ecopy(e);
        }
        else
        {
            return ecopy(e);
        }
    }
    else
    {
        return ecopy(e);
    }
    return ecopy(e);
}
