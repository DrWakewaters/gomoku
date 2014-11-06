#ifndef GAME_H
#define GAME_H

#include <iostream>
#include <vector>
#include "Player.h"
#include "AI.h"

class Game {
public:
	Game(int maxSearchDepth, int boardHeight, int boardWidth, double maxSearchTime, bool blackIsHuman, bool whiteIsHuman, bool printAIInformation);
	~Game();
	bool makeBlackPly();
	bool makeWhitePly();
	Ply findStartingAIPly();
	bool boardIsFull();
	int getNumberOfPliesMade();
	void printBoard();
private:
	bool makeHumanPly(bool isBlack);
	bool makeAIPly(bool isBlack);
	Board *board;
	Player *blackPlayer;
	Player *whitePlayer;
	int8_t maxSearchDepth;
	int numberOfPliesMade;
	double maxSearchTime;
	bool printAIInformation;
};

#endif
