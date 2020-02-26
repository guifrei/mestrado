\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{cht20}
\LoadClass[a4paper]{article}
\RequirePackage[margin=2cm]{geometry}
\RequirePackage{titlesec}
\RequirePackage{parskip}
\RequirePackage{times}
\RequirePackage{graphicx}
\RequirePackage{authblk}
\RequirePackage[square]{natbib}

\def\pnumber#1{\def\@pnumber{#1}}
\def\ptitle#1{\def\@ptitle{#1}}
\def\pauthor#1{\def\@pauthor{#1}}
\def\pabstract#1{\def\@pabstract{#1}}

\def\maketitle{
  {\flushright
    \fontfamily{arial}\selectfont
    \fontsize{10pt}{11pt}\selectfont\sffamily
    Proceedings of CHT-20\\
    ICHMT International Symposium on Advances in Computational Heat Transfer\\ \ \\
    August 16 - 20, 2020, Rio de Janeiro, Brazil\\ \ \\
    \@pnumber\par
  }
  \par\null\par\null\par
  {\centering
    \textbf{\@ptitle}
    \par\null\par
    \@author\par
  }
  \par\null\par\null\par\null\par\null\par
}

\def\makeabstract{
  \textbf{ABSTRACT}~~\@pabstract
}

\titleformat{\section}
  {\normalsize\centering\bfseries\MakeUppercase}
  {}{0em}
  {}
  
\titleformat{\subsection}[runin]
  {\normalsize\bfseries}
  {}{0em}
  {}

\titleformat{\subsubsection}[runin]
  {\normalsize\itshape}
  {}{0em}
  {\underline}
  
\def\@fnsymbol#1{\ensuremath{\ifcase#1\or *\or \dagger\or \ddagger\or
   \mathsection\or \mathparagraph\or \|\or **\or \dagger\dagger
   \or \ddagger\ddagger \else\@ctrerr\fi}}
   
\renewcommand{\thefootnote}{\fnsymbol{footnote}}
\renewcommand{\footnotesize}{\normalsize}

