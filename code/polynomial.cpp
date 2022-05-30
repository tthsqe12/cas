#include "globalstate.h"
#include "sudcode.h"
#include "polynomial.h"
#include "ex_print.h"
#include "code.h"


void eval_poly_ex(std::vector<wex> & stack, uint8_t * prog, size_t progsize)
{
	size_t i = 0;
	while (i < progsize)
	{
//std::cout << " i = " << i << " progsize = " << progsize << std::endl;
		uint8_t t = prog[i++];
		size_t a, b, c;
		switch (t)
		{
			case COMPILEDPOLY_MUL:
				a = prog[i++];
				b = prog[i++];
				c = prog[i++];
//std::cout << "&"<<a << " = " << "&"<<b << " * " << "&"<<c << std::endl;
				stack[a].reset(ex_mulr(stack[b].get(), stack[c].get()));
//std::cout << "&"<<a << ": " << ex_tostring(stack[a].get()) << std::endl;
				break;
			case COMPILEDPOLY_ADD:
				a = prog[i++];
				b = prog[i++];
				c = prog[i++];
//std::cout << "&"<<a << " = " << "&"<<b << " + " << "&"<<c << std::endl;
				stack[a].reset(ex_addr(stack[b].get(), stack[c].get()));
//std::cout << "&"<<a << ": " << ex_tostring(stack[a].get()) << std::endl;
				break;
			case COMPILEDPOLY_IMM:
			{
				a = prog[i++];
				ulong n = 0;
				for (int j = 0; j < 8; j++)
				{
					n += ((ulong)(prog[i++])) << (8*j);
				}
//std::cout << "&"<<a << " = " << slong(n) << std::endl;
				stack[a].reset(emake_int_si(slong(n)));
//std::cout << "&"<<a << ": " << ex_tostring(stack[a].get()) << std::endl;
				break;
			}
			default:
				assert(false);
		}
	}
}


void eval_poly_fmpq(std::vector<xfmpq_t> & stack, uint8_t * prog, size_t progsize)
{
	size_t i = 0;
	while (i < progsize)
	{
//std::cout << " i = " << i << " progsize = " << progsize << std::endl;
		uint8_t t = prog[i++];
		size_t a, b, c;
		switch (t)
		{
			case COMPILEDPOLY_MUL:
				a = prog[i++];
				b = prog[i++];
				c = prog[i++];
//std::cout << "&"<<a << " = " << "&"<<b << " * " << "&"<<c << std::endl;
				fmpq_mul(stack[a].data, stack[b].data, stack[c].data);
//std::cout << "&"<<a << ": " << stack[a].tostring() << std::endl;
				break;
			case COMPILEDPOLY_ADD:
				a = prog[i++];
				b = prog[i++];
				c = prog[i++];
//std::cout << "&"<<a << " = " << "&"<<b << " + " << "&"<<c << std::endl;
				fmpq_add(stack[a].data, stack[b].data, stack[c].data);
//std::cout << "&"<<a << ": " << stack[a].tostring() << std::endl;
				break;
			case COMPILEDPOLY_IMM:
			{
				a = prog[i++];
				ulong n = 0;
				for (int j = 0; j < 8; j++)
				{
					n += ((ulong)(prog[i++])) << (8*j);
				}
//std::cout << "&"<<a << " = " << slong(n) << std::endl;
				fmpq_set_si(stack[a].data, slong(n), 1);
//std::cout << "&"<<a << ": " << stack[a].tostring() << std::endl;
				break;
			}
			default:
				assert(false);
		}
	}
}



void compile_poly(compiledpoly & prog, const poly & p)
{
	size_t nvars = p.vars.size();
	prog.clear();
	prog.nvars = nvars;
	prog.append_imm(nvars, 0);
	for (size_t i = 0; i < p.size(); i++)
	{
		assert(eis_int(p.coeffs[i].get()));
		assert(fmpz_fits_si(eint_data(p.coeffs[i].get())));
		prog.append_imm(nvars + 1, fmpz_get_si(p.coeffs[i].int_data()));
		for (size_t j = 0; j < nvars; j++)
		{
			assert(fmpz_fits_si(p.exps.data()[i*nvars + j].data));
			slong n = fmpz_get_si(p.exps.data()[i*nvars + j].data);
			while (n-- > 0)
			{
				prog.append_mul(nvars + 1, nvars +1, j);
			}
		}
		prog.append_add(nvars, nvars, nvars + 1);
	}
}


int monomial_cmp(const xfmpz_t * A, const xfmpz_t * B, size_t nvars)
{
    for (size_t i = 0; i < nvars; i++)
	{
        int cmp = fmpz_cmp(A[i].data, B[i].data);
        if (cmp != 0)
            return cmp;
    }
    return 0;
}

int monomial_gt(const xfmpz_t * A, const xfmpz_t * B, size_t nvars)
{
    for (size_t i = 0; i < nvars; i++)
    {
        int cmp = fmpz_cmp(A[i].data, B[i].data);
        if (cmp != 0)
            return cmp > 0;
    }
    return false;
}

bool monomial_equal(const xfmpz_t * A, const xfmpz_t * B, size_t nvars)
{
    for (size_t i = 0; i < nvars; i++)
	{
        if (!fmpz_equal(A[i].data, B[i].data))
			return false;
    }
    return true;
}


void monomial_add(xfmpz_t * A, const xfmpz_t * B, const xfmpz_t * C, size_t nvars)
{
	for (size_t i = 0; i < nvars; i++)
	{
		fmpz_add(A[i].data, B[i].data, C[i].data);
	}
}

void monomial_swap(xfmpz_t * A, xfmpz_t * B, size_t nvars)
{
	for (size_t i = 0; i < nvars; i++)
	{
		fmpz_swap(A[i].data, B[i].data);
	}
}

ex poly::get_coeff(const xfmpz_t * e)
{
    size_t nvars = vars.size();
    for (size_t i = 0; i < coeffs.size(); i++)
        if (monomial_equal(exps.data() + i*nvars, e, nvars))
            return coeffs[i].copy();
    return emake_cint(0);
}



class poly_term_cmp {
public:
    size_t nvars;
    const xfmpz_t * expdata;

    poly_term_cmp(const poly & p)
    {
        nvars = p.vars.size();
        expdata = p.exps.data();
    }

    bool operator()(size_t a, size_t b) const
    {
        return monomial_gt(expdata + nvars*a, expdata + nvars*b, nvars);
    }   
};

void poly_sort_terms(poly & A)
{
	size_t nvars = A.vars.size();

    if (A.coeffs.size() < 2)
        return;
/*
std::cout << "starting sorting" << std::endl;
std::cout << "A: " << A.tostring() << std::endl;
*/
    std::vector<size_t> ord(A.coeffs.size(), 0);

    for (size_t i = 0; i < A.coeffs.size(); i++)
        ord[i] = i;
/*
std::cout << "filled in ord" << std::endl;
    for (size_t i = 0; i < A.coeffs.size(); i++)
        std::cout << ord[i] << std::endl;
*/

    poly_term_cmp comparer(A);
    std::sort(ord.begin(), ord.end(), comparer);
/*
std::cout << "sorted ord" << std::endl;
    for (size_t i = 0; i < A.coeffs.size(); i++)
        std::cout << ord[i] << std::endl;
*/
    poly B(A.vars);
    for (size_t i = 0; i < A.coeffs.size(); i++)
        B.pushtermx(A.coeffs[ord[i]].copy(), A.exps.data() + nvars*ord[i]);
    A.swap(B);
/*
std::cout << "done sorting" << std::endl;

std::cout << "A: " << A.tostring() << std::endl;

	bool inorder = false;
	while (!inorder)
	{
		inorder = true;
		for (size_t i = 0; i + 1 < A.size(); i++)
		{
			if (monomial_cmp(A.exps.data() + i*nvars, A.exps.data() + (i+1)*nvars, nvars) < 0)
			{
				monomial_swap(A.exps.data() + i*nvars, A.exps.data() + (i+1)*nvars, nvars);
				A.coeffs[i].swap(A.coeffs[i+1]);
				inorder = false;
                assert(false);
			}
		}
	}
*/
}

void poly_combine_like_terms(poly & A)
{
	size_t nvars = A.vars.size();
	slong out = -1;
	for (slong in = 0; in < A.size(); in++)
	{
		assert(out < in);
		if (out > 0 && monomial_equal(A.exps.data() + out*nvars, A.exps.data() + in*nvars, nvars))
		{
			A.coeffs[out].reset(ex_addr(A.coeffs[out].get(), A.coeffs[in].get()));
		}
		else
		{
			if (out < 0 || !eis_zero(A.coeffs[out].get()))
			{
				out++;
			}
			if (out != in)
			{
				monomial_swap(A.exps.data() + out*nvars, A.exps.data() + in*nvars, nvars);
				A.coeffs[out].swap(A.coeffs[in]);
			}
		}
	}
	if (out < 0 || !eis_zero(A.coeffs[out].get()))
	{
		out++;
	}
	A.resize(out);
}

void poly_add_vars_match(poly & a, const poly & b, const poly & c)
{
	assert(&a != &b);
	assert(&a != &c);
	size_t nvars = b.vars.size();
	assert(a.vars.size() == nvars);
	assert(c.vars.size() == nvars);
	uex t;

	size_t i = 0, j = 0;
	a.clear();
	while (i < b.size() && j < c.size())
	{
		int cmp = monomial_cmp(b.exps.data() + i*nvars, c.exps.data() + j*nvars, nvars);
		if (cmp > 0)
		{
			a.pushtermr(b.coeffs[i].get(), b.exps.data() + i*nvars);
			i++;
		}
		else if (cmp < 0)
		{
			a.pushtermr(c.coeffs[j].get(), c.exps.data() + j*nvars);
			j++;
		}
		else
		{
			t.setz(ex_addr(b.coeffs[i].get(), c.coeffs[j].get()));
			if (eis_zero(t.get()))
				eclear(t.release());
			else
				a.pushtermx(t.release(), b.exps.data() + i*nvars);
			i++;
			j++;
		}
	}
	while (i < b.size())
	{
		a.pushtermr(b.coeffs[i].get(), b.exps.data() + i*nvars);
		i++;
	}
	while (j < c.size())
	{
		a.pushtermr(c.coeffs[j].get(), c.exps.data() + j*nvars);
		j++;
	}
}

void poly_sub_vars_match(poly & a, const poly & b, const poly & c)
{
	assert(&a != &b);
	assert(&a != &c);
	size_t nvars = b.vars.size();
	assert(a.vars.size() == nvars);
	assert(c.vars.size() == nvars);
	uex t;

	size_t i = 0, j = 0;
	a.clear();
	while (i < b.size() && j < c.size())
	{
		int cmp = monomial_cmp(b.exps.data() + i*nvars, c.exps.data() + j*nvars, nvars);
		if (cmp > 0)
		{
			a.pushtermr(b.coeffs[i].get(), b.exps.data() + i*nvars);
			i++;
		}
		else if (cmp < 0)
		{
			a.pushtermx(ex_negate(c.coeffs[j].copy()), c.exps.data() + j*nvars);
			j++;
		}
		else
		{
			t.setz(ex_subr(b.coeffs[i].get(), c.coeffs[j].get()));
			if (eis_zero(t.get()))
				eclear(t.release());
			else
				a.pushtermx(t.release(), b.exps.data() + i*nvars);
			i++;
			j++;
		}
	}
	while (i < b.size())
	{
		a.pushtermr(b.coeffs[i].get(), b.exps.data() + i*nvars);
		i++;
	}
	while (j < c.size())
	{
		a.pushtermx(ex_negate(c.coeffs[j].copy()), c.exps.data() + j*nvars);
		j++;
	}
}

void poly_mul_vars_match(poly & a, const poly & b, const poly & c)
{
	assert(&a != &b);
	assert(&a != &c);
	size_t nvars = b.vars.size();
	assert(a.vars.size() == nvars);
	assert(c.vars.size() == nvars);
	uex t;
	std::vector<xfmpz_t> te;
	te.resize(nvars);
//std::cout << "poly mul b size: " << b.size() << std::endl;
//std::cout << "poly mul c size: " << c.size() << std::endl;
	a.clear();
	for (size_t i = 0; i < b.size(); i++)
	{
		for (size_t j = 0; j < c.size(); j++)
		{
			t.setz(ex_mulr(b.coeffs[i].get(), c.coeffs[j].get()));
			monomial_add(te.data(), b.exps.data() + i*nvars, c.exps.data() + j*nvars, nvars);
			a.pushtermx(t.release(), te.data());
		}
	}
//std::cout << "poly mul a size: " << a.size() << std::endl;

	poly_sort_terms(a);
//std::cout << "poly mul a size: " << a.size() << std::endl;
	poly_combine_like_terms(a);
//std::cout << "poly mul a size: " << a.size() << std::endl;
}

int poly_pow(poly & a, const poly & b, ulong power)
{
	a.vars = b.vars;
	a.constant(eget_cint(1));
	poly t(0);
	t.vars = b.vars;
	for (ulong i=0; i < power; i++)
	{
		poly_mul_vars_match(t, a, b);
		a.swap(t);
	}
	return 0;
}


/* check if e is free of stuff comparing equal to elements of vlist */
static bool _freeq(er e, er vlist)
{
	if (!eis_node(e))
	{
		for (size_t i = 1; i <= elength(vlist); i++)
		{
			if (ex_same(e, echild(vlist,i)))
				return false;
		}
		return true;
	}
	else
	{
		for (size_t i = 0; i <= elength(e); i++)
		{
			if (!_freeq(echild(e,i), vlist))
				return false;
		}
		return true;
	}

}

bool ex_to_polynomial(poly & p, er expr, er vlist)
{
//std::cout << "ex_to_polynomial called" << std::endl;
//std::cout << "expr: " << ex_tostring(expr) << std::endl;
//std::cout << "vlist: " << ex_tostring(vlist) << std::endl;

	size_t nvars = p.vars.size();
	assert(p.vars.size() == elength(vlist));
	if (!eis_node(expr))
	{
		for (size_t i = 0; i < nvars; i++)
		{
			if (ex_same(expr, echild(vlist,i+1)))
			{
				p.generator(i);
				return true;
			}
		}
		p.constant(expr);
		return true;
	}
	else
	{
		poly q(nvars);
		poly t(nvars);
		if (ehas_head_sym(expr, gs.sym_sPlus.get()))
		{
			p.clear();
			for (size_t i = 0; i < elength(expr); i++)
			{
				if (!ex_to_polynomial(q, echild(expr,i+1), vlist))
				{
					return false;
				}
				poly_add_vars_match(t, p, q);
				p.swap(t);
			}
			return true;
		}
		else if (ehas_head_sym(expr, gs.sym_sMinus.get()))
		{
			p.clear();
			for (size_t i = 1; i <= elength(expr); i++)
			{
				if (!ex_to_polynomial(q, echild(expr,i), vlist))
				{
					return false;
				}
                if (i == 1)
                    poly_add_vars_match(t, p, q);
                else
                    poly_sub_vars_match(t, p, q);
				p.swap(t);
			}
			return true;
		}
		else if (ehas_head_sym(expr, gs.sym_sTimes.get()))
		{
			p.constant(eget_cint(1));
			for (size_t i = 1; i <= elength(expr); i++)
			{
				if (!ex_to_polynomial(q, echild(expr,i), vlist))
				{
					return false;
				}
				poly_mul_vars_match(t, p, q);
				p.swap(t);
			}
			return true;
		}
		else if (ehas_head_sym_length(expr, gs.sym_sPower.get(), 2)
					&& eis_int(echild(expr,2))
					&& fmpz_sgn(eint_data(echild(expr,2))) >= 0)
		{
			if (!ex_to_polynomial(q, echild(expr,1), vlist))
			{
				return false;
			}
			assert(fmpz_abs_fits_ui(eint_data(echild(expr,2))));
			poly_pow(p, q, fmpz_get_ui(eint_data(echild(expr,2))));
			return true;
		}
		else
		{
            for (size_t i = 1; i <= nvars; i++)
            {
                if (ex_same(expr, echild(vlist,i)))
                {
                    p.generator(i-1);
                    return true;
                }
            }
			if (_freeq(expr, vlist))
			{
				p.constant(expr);
				return true;
			}
			else
			{
				return false;
			}
		}
    }
}
