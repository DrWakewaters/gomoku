#ifndef GAME_H
#define GAME_H

#include <chrono>
#include <iostream>
#include <random>
#include <vector>
#include "AI.h"
#include "Player.h"

class Game {
public:
	Game(int maxSearchDepth, int boardHeight, int boardWidth, double maxSearchTime, bool blackIsHuman, bool whiteIsHuman, bool printAIInformation);
	~Game();
	bool makeBlackPly();
	bool makeWhitePly();
	Position findStartingAIPly();
	bool boardIsFull();
	int getNumberOfPliesMade();
	void printBoard();
	Board *board;
private:
	bool makeHumanPly(bool isBlack);
	bool makeAIPly(bool isBlack);
	Player *blackPlayer;
	Player *whitePlayer;
	int8_t maxSearchDepth;
	int numberOfPliesMade;
	double maxSearchTime;
	bool printAIInformation;
};

#endif
