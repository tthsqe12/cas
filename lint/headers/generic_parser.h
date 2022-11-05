#pragma once

#define PREC_LOWEST    0
#define PREC_PLUS      1
#define PREC_MINUS     1
#define PREC_TIMES     2
#define PREC_DIVIDES   2
#define PREC_UPLUS     3
#define PREC_UMINUS    3
#define PREC_POWER     4
#define PREC_HIGHEST 255

#define OP_TIMES    0
#define OP_PLUS     1
#define OP_MINUS    2
#define OP_DIVIDES  3
#define OP_LROUND   4
#define FIX_INFIX       0
#define FIX_PREFIX      1
#define FIX_POSTFIX     2
#define FIX_MATCHFIX    3

#define _is_op(a) ((a) >= 0)
#define _op_make(name, fix, prec) (((prec) << 10) + ((fix) << 8) + ((name) << 0))
#define _op_prec(a)  ((slimb)( (ulong)(a) >> 10      ))
#define _op_fix(a)   ((slimb)(((ulong)(a) >> 8) & 3  ))
#define _op_name(a)  ((slimb)(((ulong)(a) >> 0) & 255))

template <typename Ring_T, typename Elem_T>
struct elem_parser {
    Ring_T* ring;
    Elem_T* res;
    std::vector<slimb> stack;
    std::vector<Elem_T> estore;
    std::vector<std::string> terminal_strings;
    std::vector<Elem_T> terminal_values;

    elem_parser(Ring_T& r) : ring(&r), res(nullptr) {}

    bool top_is_expr() {return !stack.empty() && !_is_op(stack.back());}

    Elem_T& top_expr() {return estore.at(-1-stack.back());}

    void push_op(slimb op)
    {
        FLINT_ASSERT(_is_op(op));
        stack.push_back(op);
    }

    /* if the top is not an expr, push the *res, otherwise fail */
    void push_expr() {
        if (top_is_expr())
            throw "syntax_error";
        stack.push_back(-1-slimb(estore.size()));
        estore.push_back(std::move(*res));
    }

    /* if the top is an expr, pop it, otherwise fail */
    void pop_expr() {
        if (!top_is_expr())
            throw "syntax_error";
        swap(*ring, *res, estore.at(-1-stack.back()));
        estore.pop_back();
        stack.pop_back();
    }

    /* if the top is an operation op, pop it, otherwise fail */
    void pop_op(slimb op)
    {
        slimb n = slimb(stack.size()) - 1;

        if (n < 0 || !_is_op(stack.at(n)) || _op_name(stack.at(n)) != op)
            throw "syntax_error";

        stack.pop_back();
    }

    /* pop ops with precedence > prec */
    void pop_prec(slimb prec)
    {
        slimb n, n1, n2, n3, p, l1, l3;

        if (stack.empty())
            throw "syntax_error";

    again:

        n = stack.size();
        if (n < 2)
            return;

        n1 = stack.at(n-1);
        n2 = stack.at(n-2);
        if (_is_op(n1) || !_is_op(n2))
            return;

        n1 = -1-n1;

        p = _op_prec(n2);
        if (p < prec)
            return;

        if (_op_fix(n2) == FIX_INFIX)
        {
            n3 = stack.at(n-3);
            FLINT_ASSERT(!_is_op(n3));
            n3 = -1 - n3;
            FLINT_ASSERT(n1 == n3 + 1);

            if (_op_name(n2) == OP_TIMES)
                mul(*ring, *res, estore.at(n3), estore.at(n1));
            else if (_op_name(n2) == OP_PLUS)
                add(*ring, *res, estore.at(n3), estore.at(n1));
            else if (_op_name(n2) == OP_MINUS)
                sub(*ring, *res, estore.at(n3), estore.at(n1));
            else if (_op_name(n2) == OP_DIVIDES)
                divexact(*ring, *res, estore.at(n3), estore.at(n1));
            else
                FLINT_ASSERT_ALWAYS(false && "_pop_stack: internal error");

            swap(*ring, estore.at(n3), *res);
            estore.pop_back();
            stack.pop_back(); stack.pop_back();

            goto again;
        }
        else if (_op_fix(n2) == FIX_PREFIX)
        {
            if (_op_name(n2) == OP_MINUS)
                neg(*ring, estore.at(n1));

            stack.at(n-2) = -1-n1;
            stack.pop_back();
            goto again;
        }
        else
        {
            return;
        }
    }

    void parse(Elem_T& ans, const char* s, slimb slen)
    {
        const char * send = s + slen;
        fmpz c;

        res = &ans;

        while (s < send)
        {
            if ('0' <= *s && *s <= '9')
            {
                const char* end = s + 1;
                while (end < send && '0' <= *end && *end <= '9')
                    end++;
                c.set_strn(s, end - s);
                s = end;

                set(*ring, *res, c);
                push_expr();
            }
            else if (*s == '^')
            {
                if (++s >= send || !('0' <= *s && *s <= '9'))
                    throw "syntax_error: bad exponent";

                const char* end = s + 1;
                while (end < send && '0' <= *end && *end <= '9')
                    end++;
                c.set_strn(s, end - s);
                s = end;
                
                pop_prec(PREC_POWER);

                if (!top_is_expr())
                    throw "syntax error";

                pow(*ring, *res, top_expr(), c);
                swap(*ring, top_expr(), *res);
            }
            else if (*s == '*')
            {
                if (!top_is_expr())
                    throw "syntax error: unary *";

                pop_prec(PREC_TIMES);
                push_op(_op_make(OP_TIMES, FIX_INFIX, PREC_TIMES));
                s++;
            }
            else if (*s == '+')
            {
                if (!top_is_expr())
                {
                    push_op(_op_make(OP_PLUS, FIX_PREFIX, PREC_UPLUS));
                }
                else
                {
                    pop_prec(PREC_PLUS);
                    push_op(_op_make(OP_PLUS, FIX_INFIX, PREC_PLUS));
                }
                s++;
            }
            else if (*s == '-')
            {
                if (!top_is_expr())
                {
                    push_op(_op_make(OP_MINUS, FIX_PREFIX, PREC_UMINUS));
                }
                else
                {
                    pop_prec(PREC_MINUS);
                    push_op(_op_make(OP_MINUS, FIX_INFIX, PREC_MINUS));
                }
                s++;
            }
            else if (*s == '/')
            {
                if (!top_is_expr())
                    throw "syntax error: unary /";

                pop_prec(PREC_DIVIDES);
                push_op(_op_make(OP_DIVIDES, FIX_INFIX, PREC_DIVIDES));
                s++;
            }
            else if (*s == ' ')
            {
                s++;
            }
            else if (*s == '(')
            {
                if (top_is_expr())
                    throw "syntax error: (";

                push_op(_op_make(OP_LROUND, FIX_MATCHFIX, PREC_LOWEST));
                s++;
            }
            else if (*s == ')')
            {
                pop_prec(PREC_LOWEST);
                pop_expr();
                pop_op(OP_LROUND);
                push_expr();
                s++;
            }
            else
            {
                for (ulimb k = 0; k < terminal_values.size(); k++)
                {
                    ulimb l = terminal_strings[k].size();
                    if (0 == strncmp(s, terminal_strings[k].c_str(), l))
                    {
                        set(*ring, *res, terminal_values[k]);
                        push_expr();
                        s += l;
                        goto continue_outer;
                    }
                }
                throw "syntax error: unknown token";
            }
        continue_outer:;
        }

        pop_prec(PREC_LOWEST);
        pop_expr();

        if (!stack.empty())
            throw "syntax error";
    }

    // the token (s, sn) have value v
    void add_terminal(typename Elem_T::source_t v, const char* s, ulimb sn)
    {
        terminal_strings.push_back(std::string(s, sn));
        terminal_values.push_back(v);
        // TODO sort by length
    }

    // for c strings
    void parse(Elem_T& x, const char* s) {parse(x, s, strlen(s));}
    void add_terminal(typename Elem_T::source_t v, const char* s) {add_terminal(v, s, strlen(s));}
};




